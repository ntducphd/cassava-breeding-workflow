# Cassava simulation methods

The implemented configuration is `configs/pilot.json`. Eight independently seeded founder replicates each evaluate cycles 0 and 1, with one intervening parental advance. The three strategies share initial populations and observation seeds within each replicate. All genetic maps, trait relationships, environmental gradients and treatment costs are assumptions.

## Genetic population

`src/uavbreed/genetics.R` uses AlphaSimR to generate phased diploid haplotypes, meiosis and recombination. Sixty-four founders produce 16 families with six progeny each. The 18 chromosomes each span 1.0 Morgan and have 64 segregating sites. Alleles start at probability 0.5 and refresh between positions with probability 1 − exp(−Δm/0.08). These synthetic linkage blocks are not physical assembly coordinates or an inferred population history.

Four additive–dominance traits represent yield, quality, water response and a generic resistance score. Each has 12 QTL per chromosome, additive variance 1 and mean 0; additive and dominance-degree correlation matrices are identity matrices. Dominance degree has mean 0.3 and variance 0.1. The 24 markers per chromosome are disjoint from the QTL. Total genetic values `gv1`–`gv4` and breeding values `bv1`–`bv4` are exported separately. Scaling uses the initial cohort's total-genetic-value SD throughout both cycles.

Twelve parents advance, with at most two per family. Clonal deployment retains total genetic performance; parent recycling uses an estimated additive-marker component, distinct from the exported true breeding values. No CMD2 allele, infection assay or epidemic is implemented. Genetic and breeding sources are listed in `references/Programme_References.bib` and `references/cassava_evidence.md`.

## Field and canopy observations

`programme.py` crosses clone IDs with three sites, three harvest/year indices, two management levels and two field replicates. With site s, harvest t and management m indexed from zero:

    w = 0.25 + 0.12s + 0.10t − 0.16m
    n = 0.4 + 0.2m + 0.05s
    p = 0.15 + 0.08s + 0.06t
    d = logistic(2p − 0.65g4 − 1)
    η = 28 + 5(g1 + 0.30g3w + 0.20g1n − 0.80w + 0.30n − 0.50d)
    Y = η + 5(B + ε)

Here gk is the standardised total genetic value `gvk`; B is a block effect with SD 0.12 and ε has SD 0.65. Cassava has neither the coffee-specific bearing term nor its permanent-plot effect. Fresh-root yield is in assumed t ha⁻¹. Root dry matter is generated as 33 + 3(g2 + 0.10w + εq), with εq SD 0.25. A bounded disease score is generated separately.

Height, cover, NDRE, thermal and senescence proxies combine genetic values, environment and noise. Missingness is 0.08 and base sensor-noise SD is 0.2. These traits are generated directly; no raw-image reconstruction or root-quality assay is simulated.

## Prediction and advancement

`models.py` builds additive and dominance features from training marker frequencies. Ridge regression (penalty 1) predicts yield, dry matter and the disease score, with environmental inputs and optional candidate-canopy traits. Median imputation, missingness indicators and centring are fitted to training rows. The parent model contains additive markers and environment, without dominance or canopy features. The two genomic strategies therefore share parental rankings and recycled cohorts.

The product score uses standardised predicted traits with weights (0.65, 0.20, −0.15); these are scenario weights. All four holdouts reserve the final harvest. `new_parents`, `new_site` and `new_parents_site` additionally exclude the specified parent pair, site or both from training. Candidate canopies remain observed at scoring. The field baseline uses past clone means and a training mean for unseen IDs. A five-year cycle is assumed equally across strategies; it does not measure accelerated cultivar release or realised additive gain.

## Selected-clone management

`deployment.py` selects the final-cycle top 20% under `genomic_uav`: 19 clones per replicate. Each supplies 600 training and 1,000 test fields with the same clone IDs represented in both sets. Training treatment is randomised within clone. For stress z drawn uniformly on [0,1], latent yield value vg and response rg:

    Y(0) = 28 + 5(vg + 0.20rg − z)
    δ = 0.90z + 0.15rg − 0.15
    Y(1) = Y(0) + 5δ

With probability 0.25, δ becomes −0.12. Stress and response proxies have noise SDs 0.2 and 0.4; observed outcome noise SD is 0.35 × 5. Learners use clone IDs, the proxies and one nuisance feature. Two gradient-boosting models use the randomised treatment arms. None, blanket treatment, sensor >0.5, learned policy and a truth-only oracle are compared on common test fields through Y(a) − 0.3 × 5 × a. The separate management experiment inside `programme.py` uses independently generated response values rather than selected clone IDs.

## Supplementary Table S6

Four prediction contrasts subtract genomic from genomic-plus-canopy RMSE divided by the assumed yield scale 5. Two cycles are averaged within each founder replicate. Three management contrasts compare net yield equivalents on the same held-out fields. For each set of eight paired differences d, S6 reports the mean, SD(d), SE = SD(d)/sqrt(8), and mean(d) ± t(0.975,7) × SE. These pointwise 95% Monte Carlo intervals are not adjusted for multiple comparisons and do not quantify field-calibration uncertainty.

The standalone calculation checks frozen numerical-source hashes and equality of parent scores and recycled genomic cohorts. Supplied outputs comprise seven contrasts and 56 replicate values. See `README.md` for installation and replay commands.
