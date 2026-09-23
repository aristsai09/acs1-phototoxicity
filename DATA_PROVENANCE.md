# Data provenance

## Experimental endpoint

The primary phototoxicity endpoint originates from the published compilation by Schmidt et al.:

Schmidt, F.; Wenzel, J.; Halland, N.; Gussregen, S.; Delafoy, L.; Czich, A. Computational Investigation of Drug Phototoxicity: Photosafety Assessment, Photo-Toxophore Identification, and Machine Learning. *Chem. Res. Toxicol.* **2019**, *32*, 2338-2352. DOI: 10.1021/acs.chemrestox.9b00338.

The final locked modeling cohort contains 450 compounds, 192 positive and 258 negative. The public reproducibility package emphasizes derived/computed materials: exact out-of-fold predictions, split/fold identifiers, state enumerations, quantum summaries, bootstraps, calibration tables, and Q-Chem results.

## Computed state data

`data/state_tables/002_02_microstates.csv` contains the enumerated ACS1 molecular-state table used in the state audit. `data/state_tables/005_04_stda_state_summary.csv` contains the completed state-level sTDA/xTB screening summaries used for trusted-state aggregation and state-variability analysis.

Automated state enumeration can generate chemically implausible forms. The manuscript therefore distinguishes all enumerated states from the trusted-state subset used in the final state-aware summaries.

## Frozen machine-learning outputs

The authoritative final OOF matrix is:

`data/frozen_predictions/17_core_700tree_repeated_scaffold_oof.csv`

It contains the exact repeat/fold identifiers and OOF probabilities for the final 10-repeat, five-fold scaffold-grouped robustness analysis. Compound-mean probabilities are distributed separately in:

`data/frozen_predictions/19_core_700tree_compound_mean_oof_predictions.csv`

## Higher-level quantum data

`data/quantum/15_qchem_validation_statistics.csv` stores the frozen validation statistics for the 24-compound higher-level panel. `data/quantum/matched_state_qchem_results.csv` contains the numerical results for the four matched molecular-state sentinel pairs. Completed Q-Chem input/output files for those matched pairs are included under `qchem/matched_state_jobs/`.

## Not mixed into the final effect-size freeze

The 99 label-permutation campaign preceded the final 700-tree effect-size freeze. Its files are preserved as a separate no-signal robustness check and are not numerically substituted for the final paired-bootstrap effect sizes.
