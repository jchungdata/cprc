# cprc 0.1.0

* Initial CRAN release of `cprc`.
* Added supervised streaming regression (`cprc_reg()`) and binary classification (`cprc_class()`) combiners with ADMM and simplex projection (`admm_cprc()`, `prox_simplex()`).
* Added K-fold cross-validation routines (`cv_cprc_reg()`, `cv_cprc_class()`, `cv_cprc_qp()`) supporting both `"1se"` and `"min"` selection rules.
* Added unsupervised constrained ridge optimal scoring (`cprc_unsupervised()`) and parametric bootstrap weighted *t*-tests (`cprc_boot_test()`).
* Added concept drift and outlier diagnostics (`concept_drift_angle()`, `extract_clusters()`, `plot_elbow()`) and synthetic stream generators (`make_seq_reg()`, `make_seq_plane()`, `make_seq_sea()`).
* Bundled benchmark datasets (`ccpp.csv`, `spam.RData`, `synth_spm.csv`) and introductory workflow vignette (`cprc-intro`).
