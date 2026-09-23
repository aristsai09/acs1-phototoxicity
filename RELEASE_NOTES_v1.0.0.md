# v1.0.0 - ACS Omega submission reproducibility archive

Initial frozen release supporting the manuscript **Excited-State Photophysics Improves Scaffold-Robust Phototoxicity Prediction and Reveals Compound-Specific State Effects**.

This release includes:

- recovered final model configuration and original ACS1 pipeline scripts;
- exact 10-repeat x 5-fold scaffold-grouped OOF predictions and fold identifiers;
- final compound-mean predictions and repeated-CV summaries;
- 10,000-resample paired-bootstrap summaries and controls;
- label-permutation null results;
- calibration and threshold-transition tables;
- trusted molecular-state and sTDA state tables;
- higher-level Q-Chem validation statistics;
- method-matched Q-Chem input/output records for the four sentinel state pairs;
- exact environment manifests;
- final publication figures and TOC graphic;
- file-level SHA-256 checksums and a reproduction map.

## Provenance caveat

The recovered base training configuration records five scaffold-CV repeats, whereas the final locked robustness artifact contains ten repeats. The historical standalone driver that extended the run to ten repeats was not recovered. To avoid reconstructing undocumented historical behavior, the deposited final repeat/fold identifiers and OOF probabilities are treated as authoritative; all reported downstream statistics regenerate from those immutable predictions.

## Version

`v1.0.0` is the intended ACS Omega submission snapshot. Any later revision should be released as a new version without modifying this frozen release.
