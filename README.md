# cprc

**Constrained Penalized Regression Combiners for Streaming Ensembles and Unsupervised Aggregation**

`cprc` implements simplex-constrained penalized regression combiners in pure R for:

- **Supervised streaming regression** (`cprc_reg()`) under concept drift using generalized fused lasso, pairwise/identity ridge, constrained least squares, and inverse-RMSE weighting.
- **Supervised streaming binary classification** (`cprc_class()`) using penalized Bernoulli negative log-likelihood on the unit probability simplex.
- **Unsupervised multi-subject aggregation** (`cprc_unsupervised()`, `plot_elbow()`, `cprc_boot_test()`) for combining heterogeneous high-dimensional parameter maps with cumulative weight elbow diagnostics and bootstrap weighted *t*-tests.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("jchungdata/cprc")
```

## Quick Example

```r
library(cprc)
set.seed(123)

dat <- make_seq_reg(n = 200, M = 12, p = 10, change_points = c(7), sigma = 1)
mods <- lapply(dat$seqdata[1:10], function(df) lm(y ~ ., data = df))
Fm_tr <- sapply(mods, function(m) predict(m, newdata = dat$seqdata[[11]]))

fit <- cprc_reg(y = dat$seqdata[[11]]$y, Fm = Fm_tr, penalty = "lasso")
print(fit)
plot(fit)
```

## Vignette

After installation with `build_vignettes = TRUE`, view the package workflow vignette via:

```r
vignette("cprc-intro", package = "cprc")
```

## Bug Reports & Issues

Please report bugs or feature requests at <https://github.com/jchungdata/cprc/issues>.
