"""Check deployment linkage, observed inputs and paired policy evaluation."""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
import pytest
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from uavbreed.deployment import (OBSERVED_FEATURES, evaluate_policies, fit_policy,
                                 load_products, model_inputs, run_deployment, simulate_fields)


@pytest.fixture(params=['maize', 'cassava'])
def case(request):
    crop = request.param
    cfg = {"management": {"training_fields": 120, "test_fields": 100,
                           "nonresponse_rate": .25, "treatment_cost_yield_units": .3},
           "phenotyping": {"sensor_noise_sd": .2},
           "crops": {crop: {"mean_yield": 28.0 if crop == 'cassava' else 7.5,
                            "yield_scale": 5.0 if crop == 'cassava' else 1.1}}}
    selected = pd.DataFrame({"genotype_id": ["g0", "g1", "g2", "g3"],
                             "gv1": [-.5, .2, 1., 1.5], "gv3": [-1., 0., .5, 1.2]})
    observed, truth = simulate_fields(selected, cfg, crop, 0, 809)
    return cfg, selected, observed, truth


def test_saved_products_match_score_ranking_and_compatibility():
    source = ROOT / "analysis/results/validated_pilot"
    if not (source / "config_used.json").exists():
        pytest.skip("Frozen programme output is not present")
    cfg = json.loads((source / "config_used.json").read_text())
    cycle = cfg["cycles"] - 1
    for crop in cfg["crops"]:
        selected, compatibility = load_products(source, cfg, 0, crop, [])
        scores = pd.read_csv(source / f"snapshots/r000/{crop}/genomic_uav_cycle{cycle}_scores.csv",
                             dtype={"genotype_id": str})
        ranked = scores.groupby("genotype_id").score.mean().sort_values(ascending=False)
        count = 6 if crop == "coffee" else max(1, int(len(ranked) * cfg['selection_fraction']))
        assert selected.genotype_id.tolist() == ranked.index[:count].tolist()
        if crop == "coffee":
            assert set(compatibility.recipient) == set(selected.genotype_id)
            assert len(compatibility) == 36


def test_independent_fields_share_only_selected_known_genotypes(case):
    _, selected, observed, _ = case
    train, test = observed[observed.split == "train"], observed[observed.split == "test"]
    assert not set(train.field_id) & set(test.field_id)
    assert set(train.genotype_id) == set(test.genotype_id) == set(selected.genotype_id)
    assert train.groupby("genotype_id").treatment.nunique().eq(2).all()


def test_genotype_truth_changes_corresponding_potential_yields(case):
    cfg, selected, observed, truth = case
    changed = selected.copy()
    changed.loc[changed.genotype_id == "g2", "gv1"] += 1
    crop = next(iter(cfg['crops']))
    scale = cfg['crops'][crop]['yield_scale']
    other_observed, other_truth = simulate_fields(changed, cfg, crop, 0, 809)
    affected = truth.genotype_id.eq("g2").to_numpy()
    np.testing.assert_allclose(other_truth.y0_truth - truth.y0_truth, affected * scale)
    np.testing.assert_allclose(other_truth.y1_truth - truth.y1_truth, affected * scale)
    pd.testing.assert_series_equal(observed.genotype_id, other_observed.genotype_id)
    response_changed = selected.copy()
    response_changed.loc[response_changed.genotype_id == "g2", "gv3"] += 1
    _, response_truth = simulate_fields(response_changed, cfg, crop, 0, 809)
    np.testing.assert_allclose(response_truth.y0_truth - truth.y0_truth, affected * .2 * scale)
    effect_change = (response_truth.y1_truth - response_truth.y0_truth) - (truth.y1_truth - truth.y0_truth)
    expected = affected * (~truth.nonresponse_truth.to_numpy()) * .15 * scale
    np.testing.assert_allclose(effect_change, expected, atol=1e-12)


def test_no_truth_inputs_or_test_target_leakage(case):
    cfg, selected, observed, truth = case
    crop = next(iter(cfg['crops']))
    cost = cfg['management']['treatment_cost_yield_units'] * cfg['crops'][crop]['yield_scale']
    levels = selected.genotype_id.tolist()
    changed = observed.copy()
    changed["gv1_truth"] = truth.gv1_truth * 999
    changed["gv3_truth"] = truth.gv3_truth * 999
    changed["y0_truth"] = -99999
    changed["y1_truth"] = 99999
    changed.loc[changed.split == "test", "yield_observed"] = 1e8
    np.testing.assert_equal(model_inputs(observed, levels), model_inputs(changed, levels))
    assert not any("truth" in name for name in OBSERVED_FEATURES)
    with threadpool_limits(limits=1):
        a = fit_policy(observed, levels, cost, 809)
        b = fit_policy(changed, levels, cost, 809)
    pd.testing.assert_frame_equal(a, b)


def test_reproducible_and_oracle_dominates_paired_policies(case):
    cfg, selected, observed, truth = case
    other, other_truth = simulate_fields(selected, cfg, next(iter(cfg['crops'])), 0, 809)
    pd.testing.assert_frame_equal(observed, other)
    pd.testing.assert_frame_equal(truth, other_truth)
    crop = next(iter(cfg['crops']))
    cost = cfg['management']['treatment_cost_yield_units'] * cfg['crops'][crop]['yield_scale']
    with threadpool_limits(limits=1):
        prediction = fit_policy(observed, selected.genotype_id.tolist(), cost, 809)
    metrics, fields = evaluate_policies(observed, truth, prediction, cost)
    oracle = metrics.set_index("policy").loc["oracle", "net_yield_equivalent"]
    assert metrics.net_yield_equivalent.le(oracle + 1e-12).all()
    assert metrics.fields.eq(100).all()
    assert (fields.net_yield_oracle >= fields.net_yield_learned_policy - 1e-12).all()
    assert set(fields.genotype_id) == set(selected.genotype_id)


def test_output_cannot_overlap_frozen_source(tmp_path):
    source = tmp_path / "frozen"
    source.mkdir()
    with pytest.raises(ValueError):
        run_deployment(ROOT, source, source / "new")


@pytest.fixture
def saved_cassava_products(tmp_path):
    base = tmp_path / 'snapshots/r000/cassava'
    state = base / 'genomic_uav_cycle1'
    state.mkdir(parents=True)
    ids = [f'g{i}' for i in range(8)]
    pd.DataFrame({'genotype_id': ids, 'score': range(8)}).to_csv(base / 'genomic_uav_cycle1_scores.csv', index=False)
    pd.DataFrame({'genotype_id': ids, 'gv1': range(8), 'gv3': np.zeros(8)}).to_csv(state / 'genetic_truth.csv', index=False)
    return tmp_path


def test_deployment_uses_source_selection_fraction(saved_cassava_products):
    cfg = {'cycles': 2, 'selection_fraction': .5}
    selected, _ = load_products(saved_cassava_products, cfg, 0, 'cassava', [])
    assert selected.genotype_id.tolist() == ['g7', 'g6', 'g5', 'g4']


@pytest.mark.parametrize('fraction', [0, -.1, 1.1, float('nan'), True, None, '0.5'])
def test_invalid_selection_fraction_is_rejected(saved_cassava_products, fraction):
    with pytest.raises(ValueError, match='selection_fraction'):
        load_products(saved_cassava_products, {'cycles': 2, 'selection_fraction': fraction}, 0, 'cassava', [])
