# Reproducibility

## Recorded environment

The publication replay used Windows, Python 3.13.13, R 4.4.3, AlphaSimR 2.1.0 and jsonlite 2.0.0 on a CPU. Python dependencies are pinned in `requirements.lock.txt`; the R bootstrap checks package versions. Platform and dependency changes can affect fitted models. No GPU or external genotype database is needed.

## Outputs and runtime

The eight-founder, two-cycle genetic pilot completed in 58.94 seconds on the recorded workstation. The S6 demonstration from saved outputs completed in 2.80 seconds. These are observed runtimes, not hardware-independent limits. Installation time depends on package access and whether compilation is needed.

`scripts/reproduce_s6.py` produces `Supplementary_Table_S6_Simulation_Contrasts.csv` (seven rows), `Supplementary_Table_S6_Simulation_Replicates.csv` (56 rows), a workbook and the calculation script. Both CSVs must match their reference files in `data/expected_s6/` byte for byte in the recorded environment. Frozen outputs use explicit CRLF/LF bytes; `.gitattributes` disables automatic newline conversion so source and data hashes survive cloning.

`verify_release.py` checks SHA256 records, genetic-engine/configuration hashes, relative deployment provenance and the S6 replay. `audit_workflow.py` checks 32 genetic states, 192 validation partitions, disjoint parents/genotypes in new-parent tests and 19 selected clones per founder replicate. Its stage-coverage table is a documented assessment of implementation scope.

## Interpretation

All field records, genotypes and outcomes in the saved pilot are synthetic. Founder replicates quantify Monte Carlo variation under one fixed generator. Site, year and management dimensions are crossed in the generated trial. The two evaluated cycles contain one parental advance; elapsed breeding time is assumed. The parent model is additive-only in both genomic strategies. The UAV variables are directly generated proxies, rather than raw imagery. Held-out management fields contain the same selected clone IDs as their training fields.

Stage-specific seedling, CET, PYT, AYT and regional advancement, flowering/seed set, cutting capacity, pathogen-tested multiplication and disease epidemics require further models or field records. `configs/research_200_replicates.json` and `configs/sensitivity_plan.json` are prospective configurations; their results are not part of the publication pilot. `phenotyping.early_fraction_of_season` is retained in the original configuration but is unused by the implementation.

## Data and rights

The repository distributes synthetic numerical records, original code, documentation and source metadata. Bibliographic snapshots and third-party full texts are excluded; URLs, access dates and assessed locations are retained in the source ledger. A missing snapshot does not prevent simulation or S6 reproduction. Dependencies keep their upstream licences. MIT applies to the original repository material; it does not relicense dependencies or source articles.
