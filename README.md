# Excited-State Photophysics for Scaffold-Robust Phototoxicity Prediction

Reproducibility repository for the ACS1 phototoxicity study:

**Excited-State Photophysics Improves Scaffold-Robust Phototoxicity Prediction and Reveals Compound-Specific State Effects**

Author: **Aris Tsai**  
Affiliation: Stanford Online High School, Stanford University  
ORCID: 0009-0004-8641-7647

## Study overview

This repository supports a computational phototoxicity study that separates the predictive value of excited-state information from the mechanistic value of explicit molecular-state resolution.

The locked primary cohort contains 450 compounds (192 phototoxicity-positive, 258 negative). Models were evaluated with 10 repeats of five-fold scaffold-grouped cross-validation. The final staged representations were:

- **A — Structure only:** Morgan fingerprint baseline.
- **C — Parent photophysics:** A plus 11 parent-state xTB/sTDA descriptors.
- **D — Microstate summary:** state-aware aggregation to 56 quantum/photophysical features.
- **E — Microstate + conformers:** 184-feature conformer-augmented representation.

Key findings reproduced by the frozen outputs in this repository include:

- Structure only: ROC-AUC ~0.745, PR-AUC ~0.668.
- Parent photophysics: ROC-AUC ~0.792, PR-AUC ~0.716.
- Microstate summary: ROC-AUC ~0.807, PR-AUC ~0.726.
- Microstate + conformers: ROC-AUC ~0.808, PR-AUC ~0.729.
- In 10,000 paired compound-level bootstraps, the largest robust ranking increment is structure -> parent photophysics (Delta ROC-AUC ~+0.0457).
- Higher-level Q-Chem validation preserves S1/T1 ordering strongly across the 24-compound validation panel (Pearson r ~0.946/~0.931).
- Matched-state calculations show strong state dependence for avobenzone and 2-ethylanthraquinone, but comparatively small higher-level differences for pradofloxacin and ketoprofen.

The central conclusion is deliberately restrained: **excited-state information is the dominant reproducible addition beyond structure, while explicit molecular-state effects are chemically meaningful but heterogeneous rather than universally accuracy-enhancing.**

## Repository layout

```text
code/
  analysis/                 Authoritative downstream statistical analysis
  original_pipeline/        Recovered ACS1 model/state-generation scripts
  reproduction/             Relative-path reproduction helpers
config/                     Frozen model/pipeline configuration
data/
  frozen_predictions/       Exact repeat/fold OOF predictions and summaries
  analysis/                 Bootstrap, calibration, transition, and figure tables
  quantum/                  Higher-level validation and matched-state tables
  state_tables/             Computed microstate and sTDA state tables
environment/                 Frozen Python/Conda environment records
figures/                     Final publication figures and TOC graphic
qchem/matched_state_jobs/    Method-matched Q-Chem inputs/outputs used for the four sentinel pairs
```

## Reproduce the final downstream statistics

Create an environment with Python, NumPy, pandas, SciPy, and scikit-learn, then run:

```bash
python code/reproduction/reproduce_downstream.py
```

The script reads the frozen compound-mean predictions and state tables from this repository and writes regenerated metrics, 10,000-resample paired bootstrap summaries, classification transitions, and trusted-state variability correlations to `reproduced_results/`.

The exact environment used during the study is recorded in `environment/`.

## Frozen-model provenance

Recovered ACS1 files confirm:

- Morgan fingerprint radius = 2
- fingerprint length = 2048 bits
- `ExtraTreesClassifier`
- `n_estimators = 700`
- `min_samples_leaf = 2`
- `class_weight = "balanced"`
- `max_features = "sqrt"`
- `n_jobs = 4`
- base random seed = 20260907
- scaffold grouping = Bemis-Murcko scaffold
- five-fold shuffled `StratifiedGroupKFold`

The recovered base pipeline configuration recorded five repeats; the final locked robustness pass contains ten repeats. The standalone driver used to extend the robustness analysis was not recovered. Therefore the **exact final repeat/fold identifiers and OOF probabilities distributed in `data/frozen_predictions/17_core_700tree_repeated_scaffold_oof.csv` are authoritative**. Downstream published statistics are regenerated from these locked predictions rather than attempting to infer missing historical split-generation code.

## Q-Chem provenance

`qchem/matched_state_jobs/` contains the completed input/output records used for the matched-state mechanistic panel:

- 2-ethylanthraquinone: jobs 4106/4107/4118/4119
- ketoprofen: jobs 4108/4109/4120/4121
- pradofloxacin rescue pair: jobs 4149/4150/4153/4154
- avobenzone rescue pair: jobs 4151/4152/4155/4156

The authoritative numerical summary is `data/quantum/matched_state_qchem_results.csv`.

## Source-data provenance

The original experimental phototoxicity compilation is from Schmidt et al., *Chemical Research in Toxicology* (2019), DOI: 10.1021/acs.chemrestox.9b00338. This repository focuses on the derived/computed data, locked predictions, and analysis outputs supporting the manuscript. See `DATA_PROVENANCE.md` for details.

## Reproducibility boundary

This release preserves every downstream number used in the final manuscript from the frozen predictions and computed state/Q-Chem tables. It does **not** claim that the historical 10-repeat training run can be reconstructed bit-for-bit from the recovered five-repeat base driver; the exact final OOF matrix is therefore treated as the immutable published artifact.

## Citation

See `CITATION.cff`. After Zenodo archival, cite the versioned DOI shown in the repository release and Zenodo record.

## Contact

Aris Tsai  
aristsai09@gmail.com  
ORCID: https://orcid.org/0009-0004-8641-7647
