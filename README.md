# Cassava breeding and genotype-management workflow

Version 2.0.0. Companion software for the cassava review, Section 6.4 and Supplementary Table S6.

The programme simulates diploid cassava families, separate parental and clonal selection, and treatment policies for selected clone IDs. The recorded pilot comprises eight founder replicates and two evaluated breeding cycles. All observations and outcomes are synthetic. [Methods](docs/METHODS_EN.md) define the generator, models and estimands; [operational workflow](docs/WORKFLOW_EN.md) maps the field stages and Vietnam requirements.

## Install

Tested with Python 3.13, R 4.4.3 and AlphaSimR 2.1.0 on Windows, using a CPU. The S6 demonstration and tests use Python. A new genetic simulation also requires R and the pinned R packages.

```powershell
python -m venv analysis/env/python
analysis/env/python/Scripts/python.exe -m pip install -r requirements.lock.txt
analysis/env/python/Scripts/python.exe -m pip install --no-deps -e .
```

For the genetic engine, put Rscript on PATH and run:

```powershell
Rscript --vanilla bootstrap_R.R
```

The bootstrap checks R 4.4.x and installs AlphaSimR 2.1.0 and jsonlite 2.0.0. Other operating systems require matching source-package installations and have not been validated here.

## Demonstration: reproduce Table S6

```powershell
analysis/env/python/Scripts/python.exe -I -X utf8 scripts/reproduce_s6.py --programme-root . --output analysis/reproduced_s6
```

Expected outputs: `Supplementary_Table_S6_Simulation_Contrasts.csv` (seven paired contrasts), `Supplementary_Table_S6_Simulation_Replicates.csv` (56 founder-replicate values), a workbook and the calculation script. Numerical reference files are in `data/expected_s6/`. Prediction contrasts compare normalised RMSE; policy contrasts compare net yield equivalents after the assumed treatment cost. The pointwise intervals use eight founder replicates. See [reproducibility](docs/REPRODUCIBILITY.md) for recorded timings and checks.

## Replay the full experiment

Use new output directories:

```powershell
analysis/env/python/Scripts/python.exe -I -X utf8 run.py --config configs/pilot.json --output analysis/replay
analysis/env/python/Scripts/python.exe -I -X utf8 run_deployment.py --source analysis/replay --output analysis/deployment_replay
analysis/env/python/Scripts/python.exe -I -X utf8 -m pytest -q
analysis/env/python/Scripts/python.exe -I -X utf8 verify_release.py
analysis/env/python/Scripts/python.exe -I -X utf8 audit_workflow.py
```

If needed, supply `--rscript` and `--r-lib` to `run.py`. The smoke configuration runs one founder replicate; the saved publication pilot uses eight.

## Files and interpretation

| Location | Contents |
|---|---|
| `src/uavbreed/` | Genetic engine, field generator, models and deployment |
| `configs/` | Pilot, smoke and prospective scenario configurations |
| `analysis/results/validated_pilot/` | Recorded pilot outputs and genetic-state records |
| `analysis/results/genotype_deployment/` | Selected-clone management outputs |
| `analysis/validation/` | Release and workflow checks |
| `scripts/reproduce_s6.py` | Calculation of Table S6 |
| `tests/` | Information-partition, accounting and configuration tests |
| `references/` | Source metadata and parameter provenance |

Files marked `TRUTH_ONLY` contain latent simulator values for generation or evaluation. Learners receive observed features. Management evaluates new fields of known selected clone IDs. Parent selection uses the same additive-marker model in the two genomic strategies. Stage-specific advancement, clean-stem multiplication, disease epidemics, raw-image processing and field calibration remain operational requirements described in the workflow.

## Citation and reuse

Nguyen, T. D. (2026). *Cassava breeding and genotype-management workflow*, version 2.0.0. GitHub. https://github.com/ntducphd/cassava-breeding-workflow/tree/v2.0.0

Citation metadata are in `CITATION.cff`. MIT covers the original software and documentation. Dependencies retain their own licences. Third-party articles, database exports and field photographs are not distributed. An archival DOI may be added after repository deposition; no DOI has been assigned in this release.
