"""Deploy selected genotype IDs into independent synthetic management fields."""
from io import BytesIO
from pathlib import Path
import hashlib
import json
import platform

import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import HistGradientBoostingRegressor
from threadpoolctl import threadpool_limits


OBSERVED_FEATURES = ("sensor_stress", "response_proxy", "nuisance_observed")
CASE = {
    "version": "selected_genotype_deployment_v2",
    "source_strategy": "genomic_uav",
    "coffee_candidates": 6,
    "baseline_response_coefficient": 0.20,
    "stress_treatment_coefficient": 0.90,
    "genotype_treatment_coefficient": 0.15,
    "treatment_intercept": -0.15,
    "nonresponse_effect": -0.12,
    "response_proxy_noise_sd": 0.40,
    "outcome_noise_sd": 0.35,
    "heuristic_stress_threshold": 0.50,
    "learner": {"max_iter": 60, "max_leaf_nodes": 7, "l2_regularization": 2},
}
CROP_CODES = {"coffee": 101, "maize": 202, "cassava": 303}
STATUS = "SYNTHETIC_KNOWN_SELECTED_GENOTYPES_NOT_FIELD_OBSERVATIONS"


def read_source(path, provenance, role, **case):
    path = Path(path)
    raw = path.read_bytes()
    root = Path(__file__).resolve().parents[2]
    relative = path.resolve().relative_to(root).as_posix() if path.resolve().is_relative_to(root) else path.name
    provenance.append({"role": role, "path": relative,
                       "sha256": hashlib.sha256(raw).hexdigest(), **case})
    return raw


def load_products(source, cfg, replicate, crop, provenance):
    """Select from saved product scores; retain the matching genetic state."""
    cycle = int(cfg["cycles"]) - 1
    strategy = CASE["source_strategy"]
    base = Path(source) / "snapshots" / f"r{replicate:03d}" / crop
    meta = {"replicate": replicate, "crop": crop, "source_cycle": cycle}
    score_path = base / f"{strategy}_cycle{cycle}_scores.csv"
    state = base / ("initial" if cycle == 0 else f"{strategy}_cycle{cycle}")
    scores = pd.read_csv(BytesIO(read_source(score_path, provenance, "product_scores", **meta)),
                         dtype={"genotype_id": str})
    truth = pd.read_csv(BytesIO(read_source(state / "genetic_truth.csv", provenance,
                                           "genetic_truth", **meta)), dtype={"genotype_id": str})
    if truth.genotype_id.duplicated().any() or not np.isfinite(truth[["gv1", "gv3"]]).all().all():
        raise ValueError("Invalid genetic truth: " + str(state))
    if set(scores.genotype_id) != set(truth.genotype_id):
        raise ValueError("Score/state genotype mismatch: " + str(base))
    ranking = scores.groupby("genotype_id").score.mean().sort_values(ascending=False)
    if not np.isfinite(ranking).all():
        raise ValueError("Nonfinite product score")
    fraction = cfg.get('selection_fraction')
    if (isinstance(fraction, bool) or not isinstance(fraction, (int, float))
            or not np.isfinite(fraction) or not 0 < fraction <= 1):
        raise ValueError('selection_fraction must be finite and in (0, 1]')
    count = CASE["coffee_candidates"] if crop == "coffee" else max(1, int(len(truth) * fraction))
    if count > len(ranking):
        raise ValueError("Insufficient product candidates")
    selected = truth.set_index("genotype_id").loc[ranking.index[:count]].reset_index()
    selected["selection_score"] = ranking.iloc[:count].to_numpy()
    selected["selection_rank"] = np.arange(1, count + 1)
    compatibility = None
    if crop == "coffee":
        path = base / f"{strategy}_cycle{cycle}_coffee_compatibility.csv"
        compatibility = pd.read_csv(BytesIO(read_source(path, provenance, "compatibility", **meta)),
                                     dtype={"recipient": str, "pollen_donor": str})
        ids = set(selected.genotype_id)
        if (set(compatibility.recipient) != ids or set(compatibility.pollen_donor) != ids
                or len(compatibility) != count ** 2
                or compatibility.duplicated(["recipient", "pollen_donor"]).any()):
            raise ValueError("Compatibility does not match selected coffee candidates")
        alleles = selected.set_index("genotype_id")[["S1", "S2"]]
        expected = [0.0 if row.recipient == row.pollen_donor else
                    np.mean([a not in set(alleles.loc[row.recipient]) for a in alleles.loc[row.pollen_donor]])
                    for row in compatibility.itertuples()]
        if not np.allclose(expected, compatibility.synthetic_compatible_pollen_fraction):
            raise ValueError("Saved compatibility fractions disagree with S alleles")
    return selected, compatibility


def case_seed(seed, replicate, crop):
    return int(np.random.SeedSequence([int(seed), int(replicate), CROP_CODES[crop], 1706]).generate_state(1)[0])


def balanced_ids(levels, count, rng):
    assigned = np.resize(np.asarray(levels, dtype=str), count)
    rng.shuffle(assigned)
    return assigned


def randomized_treatment(ids, rng):
    treatment = np.empty(len(ids), dtype=int)
    for ident in np.unique(ids):
        where = np.flatnonzero(ids == ident)
        if len(where) < 2:
            raise ValueError("At least two fields per genotype are needed for randomization")
        first = rng.integers(0, 2)
        allocation = np.resize(np.array([first, 1 - first]), len(where))
        rng.shuffle(allocation)
        treatment[where] = allocation
    return treatment


def simulate_fields(selected, cfg, crop, replicate, seed):
    """Generate observed measurements separately from latent potential outcomes."""
    rng = np.random.default_rng(seed)
    mgmt = cfg["management"]
    c = cfg["crops"][crop]
    levels = selected.genotype_id.astype(str).tolist()
    n_train, n_test = int(mgmt["training_fields"]), int(mgmt["test_fields"])
    if min(n_train, n_test) < 2 * len(levels):
        raise ValueError("Too few fields to represent each selected genotype in both sets")
    ids = np.r_[balanced_ids(levels, n_train, rng), balanced_ids(levels, n_test, rng)]
    split = np.r_[np.repeat("train", n_train), np.repeat("test", n_test)]
    field_ids = [f"r{replicate:03d}_{crop}_{part}_{i:05d}" for i, part in enumerate(split)]
    g = selected.set_index("genotype_id").loc[ids]
    stress = rng.uniform(0, 1, len(ids))
    response = g.gv3.to_numpy()
    y0 = c["mean_yield"] + c["yield_scale"] * (
        g.gv1.to_numpy() + CASE["baseline_response_coefficient"] * response - stress)
    effect = (CASE["stress_treatment_coefficient"] * stress
              + CASE["genotype_treatment_coefficient"] * response + CASE["treatment_intercept"])
    nonresponse = rng.random(len(ids)) < mgmt["nonresponse_rate"]
    effect[nonresponse] = CASE["nonresponse_effect"]
    y1 = y0 + c["yield_scale"] * effect
    sensor = stress + rng.normal(0, cfg["phenotyping"]["sensor_noise_sd"], len(ids))
    proxy = response + rng.normal(0, CASE["response_proxy_noise_sd"], len(ids))
    nuisance = rng.normal(0, 1, len(ids))
    treatment = np.r_[randomized_treatment(ids[:n_train], rng), randomized_treatment(ids[n_train:], rng)]
    observed_yield = np.where(treatment == 1, y1, y0) + rng.normal(
        0, c["yield_scale"] * CASE["outcome_noise_sd"], len(ids))
    observed = pd.DataFrame({"replicate": replicate, "crop": crop, "field_id": field_ids,
                             "split": split, "genotype_id": ids, "sensor_stress": sensor,
                             "response_proxy": proxy, "nuisance_observed": nuisance,
                             "treatment": treatment, "yield_observed": observed_yield})
    truth = observed[["replicate", "crop", "field_id", "split", "genotype_id"]].copy()
    truth["gv1_truth"] = g.gv1.to_numpy()
    truth["gv3_truth"] = response
    truth["stress_truth"] = stress
    truth["nonresponse_truth"] = nonresponse
    truth["y0_truth"] = y0
    truth["y1_truth"] = y1
    return observed, truth


def model_inputs(observed, levels):
    """Use only known IDs and the three declared noisy observed features."""
    ids = observed.genotype_id.astype(str).to_numpy()
    if not set(ids).issubset(set(levels)):
        raise ValueError("Deployment case does not predict new genotype IDs")
    indicator = ids[:, None] == np.asarray(levels)[None, :]
    return np.c_[indicator.astype(float), observed.loc[:, list(OBSERVED_FEATURES)].to_numpy(float)]


def fit_policy(observed, levels, cost, seed):
    train = observed.split.eq("train").to_numpy()
    test = observed.split.eq("test").to_numpy()
    if not train.any() or not test.any() or observed.field_id.duplicated().any():
        raise ValueError("Invalid independent field partition")
    if set(observed.loc[train, "genotype_id"]) != set(levels) or set(observed.loc[test, "genotype_id"]) != set(levels):
        raise ValueError("Each selected genotype must occur in training and test fields")
    x = model_inputs(observed, levels)
    predictions = []
    for action in (0, 1):
        use = train & observed.treatment.eq(action).to_numpy()
        model = HistGradientBoostingRegressor(**CASE["learner"], random_state=seed)
        model.fit(x[use], observed.loc[use, "yield_observed"].to_numpy())
        predictions.append(model.predict(x[test]))
    result = observed.loc[test, ["field_id", "genotype_id"]].copy()
    result["predicted_y0"] = predictions[0]
    result["predicted_y1"] = predictions[1]
    result["learned_action"] = predictions[1] - predictions[0] > cost
    return result


def evaluate_policies(observed, truth, prediction, cost):
    test = observed.loc[observed.split.eq("test")].merge(
        prediction, on=["field_id", "genotype_id"], validate="one_to_one")
    hidden = truth.loc[truth.split.eq("test")]
    test = test.merge(hidden.drop(columns=["replicate", "crop", "split"]),
                      on=["field_id", "genotype_id"], validate="one_to_one")
    if len(test) != len(hidden):
        raise ValueError("Missing held-out evaluation fields")
    y0, y1 = test.y0_truth.to_numpy(), test.y1_truth.to_numpy()
    decisions = {"none": np.zeros(len(test), dtype=bool), "all": np.ones(len(test), dtype=bool),
                 "uav_rule": test.sensor_stress.to_numpy() > CASE["heuristic_stress_threshold"],
                 "learned_policy": test.learned_action.to_numpy(), "oracle": y1 - y0 > cost}
    records = []
    oracle = np.maximum(y0, y1 - cost).mean()
    for name, action in decisions.items():
        gross = np.where(action, y1, y0)
        utility = gross - cost * action
        records.append({"policy": name, "fields": len(test), "genotypes": test.genotype_id.nunique(),
                        "gross_yield": float(gross.mean()), "net_yield_equivalent": float(utility.mean()),
                        "treated_fraction": float(action.mean()), "cost_per_treated_field": cost,
                        "regret_vs_oracle": float(oracle - utility.mean()),
                        "delta_vs_none": float(utility.mean() - y0.mean()), "data_status": STATUS})
        test[f"action_{name}"] = action
        test[f"net_yield_{name}"] = utility
    return pd.DataFrame(records), test


def run_deployment(root, source, output):
    root, source, output = Path(root).resolve(), Path(source).resolve(), Path(output).resolve()
    if source == output or source in output.parents or output in source.parents:
        raise ValueError("Output must be separate from the frozen source directory")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Choose an empty output directory: " + str(output))
    provenance = []
    cfg_raw = read_source(source / "config_used.json", provenance, "source_configuration")
    manifest_raw = read_source(source / "run_manifest.json", provenance, "source_run_manifest")
    cfg = json.loads(cfg_raw)
    manifest = json.loads(manifest_raw)
    if cfg["replicates"] != manifest["replicates"] or cfg["cycles"] != manifest["cycles"]:
        raise ValueError("Source configuration and manifest disagree")
    if CASE["source_strategy"] not in cfg["strategies"] or not set(cfg["crops"]).issubset(CROP_CODES):
        raise ValueError("Unsupported source strategy or crop")
    metrics, candidates, compatibilities, observations, truths, evaluations, seeds = [], [], [], [], [], [], []
    with threadpool_limits(limits=1):
        for rep in range(cfg["replicates"]):
            for crop in cfg["crops"]:
                selected, compatibility = load_products(source, cfg, rep, crop, provenance)
                seed = case_seed(cfg["seed"], rep, crop)
                seeds.append({"replicate": rep, "crop": crop, "seed": seed})
                data, truth = simulate_fields(selected, cfg, crop, rep, seed)
                cost = cfg["management"]["treatment_cost_yield_units"] * cfg["crops"][crop]["yield_scale"]
                prediction = fit_policy(data, selected.genotype_id.tolist(), cost, seed)
                summary, evaluated = evaluate_policies(data, truth, prediction, cost)
                for table in (summary, selected, evaluated):
                    table["replicate"] = rep
                    table["crop"] = crop
                    table["source_cycle"] = cfg["cycles"] - 1
                    table["source_strategy"] = CASE["source_strategy"]
                    table["yield_unit"] = cfg["crops"][crop]["yield_unit"]
                metrics.append(summary)
                candidates.append(selected)
                observations.append(data)
                truths.append(truth)
                evaluations.append(evaluated)
                if compatibility is not None:
                    compatibility["replicate"] = rep
                    compatibility["crop"] = crop
                    compatibilities.append(compatibility)
                print(f"Completed deployment replicate {rep + 1}/{cfg['replicates']}: {crop}", flush=True)
    output.mkdir(parents=True, exist_ok=True)
    tables = {"policy_metrics.csv": metrics, "selected_candidates_TRUTH_ONLY.csv": candidates,
              "observations_SYNTHETIC.csv": observations, "field_truth_TRUTH_ONLY.csv": truths,
              "test_field_policies_TRUTH_ONLY.csv": evaluations}
    if compatibilities:
        tables["coffee_selected_compatibility.csv"] = compatibilities
    for name, parts in tables.items():
        pd.concat(parts, ignore_index=True).to_csv(output / name, index=False)
    (output / "source_config_used.json").write_bytes(cfg_raw)
    (output / "source_run_manifest.json").write_bytes(manifest_raw)
    resolved = {"data_status": STATUS, "case": {**CASE, 'product_fraction': cfg['selection_fraction']}, "management": cfg["management"],
                "phenotyping": {"sensor_noise_sd": cfg["phenotyping"]["sensor_noise_sd"]},
                "crops": cfg["crops"], "seeds": seeds}
    (output / "deployment_config.json").write_text(json.dumps(resolved, indent=2) + "\n", encoding="utf-8")
    for relative in ["src/uavbreed/deployment.py", "run_deployment.py", "references/maize_sources.json",
                     "references/coffee_sources.json", "references/cassava_sources.json", "docs/DEPLOYMENT_CASE.md"]:
        path = root / relative
        if path.exists():
            read_source(path, provenance, "implementation_or_conceptual_source")
    (output / "source_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    source_label = source.relative_to(root).as_posix() if source.is_relative_to(root) else 'external_source'
    result = {"data_status": STATUS, "source_directory": source_label, "replicates": cfg["replicates"],
              "crops": list(cfg["crops"]), "source_cycle": cfg["cycles"] - 1,
              "interpretation": "Independent synthetic held-out fields of known selected genotype IDs; no new-genotype validation, measured fruit set or farm feedback.",
              "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
              "scikit_learn": sklearn.__version__,
              "outputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.iterdir()) if p.is_file()}}
    (output / "run_manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return output
