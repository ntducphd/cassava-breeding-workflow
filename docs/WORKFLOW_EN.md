# Operational cassava breeding workflow

Source assessment: 19 November 2026. The complete source-to-stage record is in `references/workflow_sources_20261119.json`; the Vietnamese account is [CASSAVA_WORKFLOW_VI.md](CASSAVA_WORKFLOW_VI.md).

Population improvement and clonal product development have separate decisions. Crossing and botanical seed generate new candidates. Parental recycling targets transmitted merit and can occur when prediction, flowering and material supply support a new crossing cycle. Clonal products proceed through seedling nursery, clonal/single-row evaluation, preliminary, advanced and uniform or regional trials, followed by farmer and processor evaluation. Stage sizes, durations and advancement thresholds depend on the product profile and target environments.

| Field component | Required records | Pilot coverage |
|---|---|---|
| Product profile and germplasm | Target market, environment, checks, accession, identity and health | Scenario traits and weights |
| Crossing and seedlings | Parents, cross ID, flowering, seed set, germination and available cuttings | Meiosis; no capacity or survival process |
| Clonal evaluation and yield trials | Clone/plot ID, design, checks, season, management, harvest age and trial quality | Balanced synthetic observation grid |
| Selection decisions | Additive parental merit; clonal yield and quality under target conditions | Separate model targets; shared genomic parent model |
| Root and health references | Harvested mass and area, dry matter, starch/quality assays and pathogen results | Directly simulated traits and generic disease score |
| Mother stock and multiplication | Tested stock, pathogen panel, sample/date, propagation batch and destination | Operational requirement |
| Farmer/processor evaluation | Field outcomes, costs, acceptance and product quality | Synthetic policy comparison of known clone IDs |
| Vietnam disclosure and rights | Applicable TCVN/TCCS, variety information and cultivation process; separate DUS route where applicable | Documentary assessment |

UAV canopy traits require independent root and laboratory references. Absence of symptoms and a negative result for a defined pathogen panel are different records. Culture-derived plants require checks appropriate to their genotype and regeneration route; CMD2-related loss of resistance has been reported in particular backgrounds.

As of the assessment date, Circular 12/2026 lists rice and maize as main crops. Resolution 66.19/2026, Annex I/A/VI, provides information self-disclosure under applicable national or institutional standards for other crops. Circular 48/2011, which issued the old cassava VCU regulation, was repealed from 5 February 2024. UPOV TG/321/1 concerns DUS, separately from agronomic or plant-health performance. Recheck the applicable legal framework at deployment, including amendments effective from January 2027.

The pilot's five-year cycle and observation-grid dimensions are scenario assumptions. `early_fraction_of_season` is retained in the historical configuration and is not an implemented flight calendar. Capacity, kinship-constrained contributions, field variance components and parameter sensitivity require calibration before operational use.
