"""Reproduce Supplementary Table S6 from the synthetic cassava pilot."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy.stats import t


CONTRASTS = "Supplementary_Table_S6_Simulation_Contrasts"
REPLICATES = "Supplementary_Table_S6_Simulation_Replicates.csv"
SCRIPT_COPY = "Supplementary_Simulation_Contrasts.py"
COHORT_COLUMNS = ("mean_genetic_value", "genetic_sd", "heterozygosity", "cohort_potential_yield")
HOLDOUTS = ("future_harvest", "new_parents", "new_site", "new_parents_site")


def sha256(blob):
    return hashlib.sha256(blob).hexdigest()


class Sources:
    """Read source files and retain their relative names and byte hashes."""

    def __init__(self, root):
        self.root = root
        self.hashes = {}

    def blob(self, relative):
        path = self.root / relative
        data = path.read_bytes()
        self.hashes[path.relative_to(self.root).as_posix()] = sha256(data)
        return data

    def json(self, relative):
        return json.loads(self.blob(relative).decode("utf-8-sig"))

    def csv(self, relative):
        return list(csv.DictReader(self.blob(relative).decode("utf-8-sig").splitlines()))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def project_directory():
    for ancestor in Path(__file__).resolve().parents:
        if (ancestor / "analysis/breeding_uav_simulation/configs/pilot.json").is_file():
            return ancestor
    return None


def programme_directory():
    for ancestor in Path(__file__).resolve().parents:
        if (ancestor / 'configs/pilot.json').is_file() and (ancestor / 'src/uavbreed/genetics.R').is_file():
            return ancestor
    project = project_directory()
    return project / 'analysis/breeding_uav_simulation' if project else None


def replicate_contrasts(cfg, prediction, policy):
    rows = []
    for crop, crop_cfg in cfg["crops"].items():
        for setting in HOLDOUTS:
            for replicate in range(cfg["replicates"]):
                means = {}
                for strategy in ("genomic_uav", "genomic"):
                    group = [r for r in prediction if r["crop"] == crop and r["validation"] == setting
                             and r["strategy"] == strategy and int(r["replicate"]) == replicate]
                    require(len(group) == cfg["cycles"] and {int(r["cycle"]) for r in group}
                            == set(range(cfg["cycles"])), "Incomplete prediction-cycle records")
                    means[strategy] = statistics.mean(float(r["rmse"]) / crop_cfg["yield_scale"] for r in group)
                rows.append({"experiment": "prediction", "crop": crop, "setting": setting,
                             "contrast": "genomic_uav_minus_genomic", "replicate": replicate,
                             "first_mean": means["genomic_uav"], "second_mean": means["genomic"],
                             "paired_difference": means["genomic_uav"] - means["genomic"],
                             "unit": "RMSE / assumed crop yield_scale (dimensionless)",
                             "assumed_yield_scale": crop_cfg["yield_scale"],
                             "source_yield_unit": crop_cfg["yield_unit"],
                             "aggregation": "cycle means within founder replicate"})
        for first, second in (("learned_policy", "all"), ("uav_rule", "all"), ("learned_policy", "uav_rule")):
            for replicate in range(cfg["replicates"]):
                group = [r for r in policy if r["crop"] == crop and int(r["replicate"]) == replicate]
                by_policy = {r["policy"]: r for r in group}
                require(len(by_policy) == len(group) == 5, "Incomplete deployment-policy records")
                require(by_policy[first]["yield_unit"] == by_policy[second]["yield_unit"]
                        == crop_cfg["yield_unit"], "Deployment and configuration units disagree")
                a, b = float(by_policy[first]["net_yield_equivalent"]), float(by_policy[second]["net_yield_equivalent"])
                rows.append({"experiment": "selected_known_genotype_deployment", "crop": crop,
                             "setting": "same_1000_held_out_fields_within_replicate",
                             "contrast": first + "_minus_" + second, "replicate": replicate,
                             "first_mean": a, "second_mean": b, "paired_difference": a-b,
                             "unit": "synthetic net yield equivalent: " + crop_cfg["yield_unit"],
                             "assumed_yield_scale": crop_cfg["yield_scale"],
                             "source_yield_unit": crop_cfg["yield_unit"],
                             "aggregation": "paired field-policy means within founder replicate"})
    return rows


def summarise(rows, count):
    keys = ("experiment", "crop", "setting", "contrast", "unit", "assumed_yield_scale", "source_yield_unit")
    groups = {}
    for row in rows:
        groups.setdefault(tuple(row[k] for k in keys), []).append(row)
    summary = []
    for key, group in groups.items():
        require(sorted(r["replicate"] for r in group) == list(range(count)), "Incomplete replicate contrast")
        differences = [r["paired_difference"] for r in group]
        mean, sd = statistics.mean(differences), statistics.stdev(differences)
        se, critical = sd / math.sqrt(count), float(t.ppf(.975, count-1))
        summary.append({**dict(zip(keys, key)), "n_replicates": count, "df": count-1,
                        "first_mean": statistics.mean(r["first_mean"] for r in group),
                        "second_mean": statistics.mean(r["second_mean"] for r in group),
                        "mean_difference": mean, "SD_MC": sd, "SE_MC": se,
                        "t_0_975_df7": critical, "ci95_low": mean-critical*se,
                        "ci95_high": mean+critical*se})
    return summary


def check_shared_parent_cohorts(sources, cfg, selection, pilot):
    cohort, parents, states = [], [], []
    for crop in cfg["crops"]:
        for replicate in range(cfg["replicates"]):
            base = pilot / "snapshots" / f"r{replicate:03d}" / crop
            for cycle in range(cfg["cycles"]):
                group = [r for r in selection if r["crop"] == crop and int(r["replicate"]) == replicate
                         and int(r["cycle"]) == cycle and r["strategy"] in ("genomic", "genomic_uav")]
                indexed = {r["strategy"]: r for r in group}
                require(len(group) == len(indexed) == 2, "Incomplete selection records")
                for column in COHORT_COLUMNS:
                    a, b = float(indexed["genomic"][column]), float(indexed["genomic_uav"][column])
                    require(a == b, "Genomic parent cohorts differ: " + column)
                    cohort.append({"crop": crop, "replicate": replicate, "cycle": cycle,
                                   "column": column, "genomic": a, "genomic_uav": b, "equal": True})
                scores = [sources.csv(base / f"{strategy}_cycle{cycle}_scores.csv")
                          for strategy in ("genomic", "genomic_uav")]
                require(len(scores[0]) == len(scores[1]) and all(
                    a["genotype_id"] == b["genotype_id"] and float(a["parent_score"]) == float(b["parent_score"])
                    for a, b in zip(*scores)), "Genomic parent scores differ")
                parents.append({"crop": crop, "replicate": replicate, "cycle": cycle,
                                "rows": len(scores[0]), "id_and_parent_score_equal": True})
                if cycle == 0:
                    continue
                for name in ("genetic_truth.csv", "markers.csv") + (("parent_markers.csv",) if crop == "maize" else ()):
                    a = sources.blob(base / f"genomic_cycle{cycle}" / name)
                    b = sources.blob(base / f"genomic_uav_cycle{cycle}" / name)
                    require(a == b, "Recycled genomic state files differ: " + name)
                    states.append({"crop": crop, "replicate": replicate, "cycle": cycle,
                                   "file": name, "byte_identical": True, "sha256": sha256(a)})
                require(sources.json(base / f"genomic_cycle{cycle}/parent_selection.json")
                        == sources.json(base / f"genomic_uav_cycle{cycle}/parent_selection.json"),
                        "Selected genomic parent IDs differ")
    return {"scope": "Shared parent scores and recycled cohorts; product scores may differ.",
            "cycle0": "Both genomic branches use the same initial state directory.",
            "selection_columns": list(COHORT_COLUMNS), "cohort_checks": cohort,
            "parent_score_checks": parents, "recycled_state_checks": states,
            "parent_selection_jsons_equal": True}


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def add_sheet(workbook, title, fields, labels, rows, note):
    sheet = workbook.create_sheet(title)
    sheet.append([note])
    sheet.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(fields))
    sheet.row_dimensions[1].height = 48
    sheet.append(labels)
    for row in rows:
        sheet.append([row[field] for field in fields])
    for cells in sheet:
        for cell in cells:
            cell.font = Font(name="Times New Roman", size=12, bold=cell.row <= 2)
            cell.alignment = Alignment(vertical="top", wrap_text=True)
            if cell.row == 2:
                cell.fill = PatternFill("solid", fgColor="DCE6F1")
            if cell.row > 2 and isinstance(cell.value, float):
                cell.number_format = "0.000000"
    sheet.freeze_panes = "D3"
    sheet.auto_filter.ref = f"A2:{get_column_letter(len(fields))}{sheet.max_row}"
    sheet.row_dimensions[2].height = 36
    for index, field in enumerate(fields, 1):
        width = 22
        if field in ("setting", "contrast", "aggregation", "experiment"):
            width = 42
        elif field in ("unit", "source_yield_unit"):
            width = 58
        elif field in ("crop", "replicate", "n_replicates"):
            width = 14
        sheet.column_dimensions[get_column_letter(index)].width = width
    for index in range(3, sheet.max_row+1):
        sheet.row_dimensions[index].height = 64
    sheet.sheet_properties.pageSetUpPr.fitToPage = True
    sheet.page_setup.orientation = "landscape"
    sheet.page_setup.paperSize = sheet.PAPERSIZE_A3
    sheet.page_setup.fitToWidth = 1
    sheet.page_setup.fitToHeight = 0
    sheet.print_title_rows = "1:2"


def write_workbook(path, summary, rows):
    workbook = Workbook()
    workbook.remove(workbook.active)
    workbook.properties.creator = ""
    workbook.properties.lastModifiedBy = ""
    workbook.properties.title = "Supplementary Table S6. Synthetic simulation contrasts"
    fields = ("crop", "setting", "contrast", "mean_difference", "SD_MC", "SE_MC", "ci95_low", "ci95_high", "n_replicates", "unit", "assumed_yield_scale")
    labels = ("Crop", "Test setting", "Paired contrast", "Mean difference", "SD across replicates", "Monte Carlo SE", "95% t lower", "95% t upper", "Replicates", "Outcome unit", "Assumed yield scale")
    add_sheet(workbook, "Summary", fields, labels, summary,
              "Synthetic pilot under one assumed data-generating process. Prediction contrasts average two cycles within replicate; all holdouts include future harvest. Intervals use eight founder replicates and t(7), not field rows.")
    raw_fields = list(rows[0])
    add_sheet(workbook, "Replicates", raw_fields, [name.replace("_", " ").capitalize() for name in raw_fields], rows,
              "One paired difference per founder replicate. Prediction: genomic-UAV minus genomic normalised RMSE (positive is worse). Deployment: net yield equivalent after the assumed treatment cost (positive is better).")
    workbook.save(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--programme-root", type=Path, help="Programme directory containing src/ and analysis/results/")
    parser.add_argument("--output", type=Path, help="Destination directory; defaults to the project supplementary directory, or the working directory")
    parser.add_argument("--audit-output", type=Path, help="Optional JSON path for source hashes and cohort-equality checks")
    args = parser.parse_args()
    project = project_directory()
    native = programme_directory()
    if args.programme_root is None and native is None:
        parser.error("--programme-root is required outside the source project")
    root = (args.programme_root or native).resolve()
    output = (args.output or (project / "analysis/results/manuscript_assets/supplementary" if project else root / 'analysis/reproduced_s6')).resolve()
    protected = [root / name for name in ['src', 'configs', 'references', 'docs', 'tests',
                 'analysis/results/validated_pilot', 'analysis/results/genotype_deployment']]
    require(root != output and all(p != output and p not in output.parents and output not in p.parents for p in protected),
            'Choose an output directory outside source files and frozen results')
    sources = Sources(root)
    pilot, deploy = Path("analysis/results/validated_pilot"), Path("analysis/results/genotype_deployment")
    cfg = sources.json(pilot / "config_used.json")
    manifest = sources.json(pilot / "run_manifest.json")
    require(cfg["replicates"] == manifest["replicates"] == 8 and cfg["cycles"] == manifest["cycles"] == 2,
            "Table S6 requires the published eight-replicate, two-cycle pilot")
    require(cfg == sources.json("configs/pilot.json"), "Frozen and programme configurations disagree")
    for relative, key in (("src/uavbreed/programme.py", "programme_sha256"), ("src/uavbreed/genetics.R", "engine_sha256"), ("configs/pilot.json", "config_sha256")):
        require(sha256(sources.blob(relative)) == manifest[key], "Frozen source hash differs: " + relative)
    for relative in ("src/uavbreed/models.py", "src/uavbreed/deployment.py", deploy / "deployment_config.json", deploy / "run_manifest.json"):
        sources.blob(relative)
    rows = replicate_contrasts(cfg, sources.csv(pilot / "prediction_metrics.csv"), sources.csv(deploy / "policy_metrics.csv"))
    summary = summarise(rows, cfg["replicates"])
    require(set(cfg["crops"]) == {"cassava"} and len(rows) == 56 and len(summary) == 7, "Unexpected Table S6 dimensions")
    equality = check_shared_parent_cohorts(sources, cfg, sources.csv(pilot / "selection_metrics.csv"), pilot)
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / (CONTRASTS + ".csv"), summary)
    write_csv(output / REPLICATES, rows)
    write_workbook(output / (CONTRASTS + ".xlsx"), summary, rows)
    script = Path(__file__).read_bytes()
    (output / SCRIPT_COPY).write_bytes(script)
    report = {"data_status": "SYNTHETIC_NOT_FIELD_OBSERVATIONS", "replicates": 8,
              "contrasts": len(summary), "replicate_rows": len(rows),
              "interval": "mean paired difference +/- t(0.975,7) * SD(paired differences)/sqrt(8)",
              "scope": "Monte Carlo uncertainty under one assumed DGP; field calibration and sensitivity uncertainty are not included.",
              "shared_parent_cohorts": equality, "input_source_sha256": sources.hashes,
              "script_sha256": sha256(script), "summary": summary,
              "output_sha256": {name: sha256((output / name).read_bytes()) for name in
                                (CONTRASTS + ".csv", REPLICATES, CONTRASTS + ".xlsx", SCRIPT_COPY)}}
    if args.audit_output:
        args.audit_output.parent.mkdir(parents=True, exist_ok=True)
        args.audit_output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"contrasts": len(summary), "replicate_rows": len(rows),
                      "source_files": len(sources.hashes), "shared_parent_cohorts_equal": True,
                      "output_files": list(report["output_sha256"])}, indent=2))


if __name__ == "__main__":
    main()
