test_that("prox_simplex computes exact Euclidean projection onto the unit simplex", {
  v <- c(-0.5, 0.2, 1.5, 0.8)
  w <- prox_simplex(v)
  # Analytical KKT solution: support {3, 4}, gamma = (1.5 + 0.8 - 1) / 2 = 0.65
  expect_equal(w, c(0, 0, 0.85, 0.15), tolerance = 1e-10)
  expect_equal(sum(w), 1, tolerance = 1e-12)
  expect_true(all(w >= 0))

  # Already on simplex
  u <- c(0.25, 0.25, 0.50)
  expect_equal(prox_simplex(u), u, tolerance = 1e-12)
})

test_that("concept_drift_angle and extract_clusters return exact diagnostic values", {
  # Uniform vector has 0-degree drift
  expect_equal(concept_drift_angle(rep(0.25, 4), unit = "degrees"), 0, tolerance = 1e-10)
  expect_equal(concept_drift_angle(rep(0.25, 4), unit = "radians"), 0, tolerance = 1e-10)

  # Single active model out of 4: cos(eta) = 0.25 / (0.5 * 1) = 0.5 -> 60 degrees
  expect_equal(concept_drift_angle(c(1, 0, 0, 0), unit = "degrees"), 60, tolerance = 1e-10)

  # Exact 2-cluster recovery
  cl <- extract_clusters(c(0.01, 0.01, 0.49, 0.49), tol = 1e-3)
  expect_equal(cl$no, 2L)
  expect_equal(sort(cl$csize), c(2L, 2L))
  expect_equal(sort(cl$unique_weights), c(0.01, 0.49), tolerance = 1e-10)
})

test_that("qp solvers satisfy simplex constraints and reduce error on known weights", {
  set.seed(101)
  n <- 100
  Fm <- cbind(rnorm(n), rnorm(n), rnorm(n))
  w_true <- c(0.0, 0.3, 0.7)
  y <- drop(Fm %*% w_true) + rnorm(n, sd = 0.05)

  sol_lse <- lse_qp(y, Fm)
  expect_equal(sum(sol_lse$w), 1, tolerance = 1e-8)
  expect_equal(sol_lse$w, w_true, tolerance = 0.05)

  sol_r1 <- ridge_qp(y, Fm, lambda = 1e6)
  expect_equal(sol_r1$w, rep(1 / 3, 3), tolerance = 1e-3)

  sol_r2 <- ridge_qp_diff(y, Fm, lambda = 1e6)
  expect_equal(sol_r2$w, rep(1 / 3, 3), tolerance = 1e-3)

  sol_wa <- weighted_avg(y, Fm)
  expect_equal(sum(sol_wa$w), 1, tolerance = 1e-10)
  expect_true(sol_wa$w[3] > sol_wa$w[1])
})

test_that("cprc_reg adapts to abrupt drift, fuses regimes, and outperforms wavg", {
  set.seed(123)
  dat <- make_seq_reg(n = 120, M = 8, p = 6, change_points = c(5), sigma = 1)
  mods <- lapply(dat$seqdata[1:6], function(df) lm(y ~ ., data = df))
  Fm_tr <- sapply(mods, function(m) predict(m, newdata = dat$seqdata[[7]]))
  Fm_te <- sapply(mods, function(m) predict(m, newdata = dat$seqdata[[8]]))
  y_tr <- dat$seqdata[[7]]$y
  y_te <- dat$seqdata[[8]]$y

  fit_lasso <- cprc_reg(y = y_tr, Fm = Fm_tr, lambda = 20, penalty = "lasso", cv = FALSE)
  fit_ridge <- cprc_reg(y = y_tr, Fm = Fm_tr, lambda = 20, penalty = "ridge", type = "pairwise", cv = FALSE)
  fit_lse   <- cprc_reg(y = y_tr, Fm = Fm_tr, penalty = "lse")
  fit_wavg  <- cprc_reg(y = y_tr, Fm = Fm_tr, penalty = "wavg")

  for (fit in list(fit_lasso, fit_ridge, fit_lse, fit_wavg)) {
    expect_s3_class(fit, "cprc_reg")
    expect_equal(sum(fit$w), 1, tolerance = 1e-6)
    expect_true(all(fit$w >= 0))
  }

  # Fused lasso should place >80% weight on post-drift models M5..M6 and beat wavg
  expect_true(sum(fit_lasso$w[5:6]) > 0.80)
  expect_true(fit_lasso$clusters$no <= 3L)
  rmse_lasso <- sqrt(mean((y_te - predict(fit_lasso, newdata = Fm_te))^2))
  rmse_wavg  <- sqrt(mean((y_te - predict(fit_wavg, newdata = Fm_te))^2))
  expect_true(rmse_lasso < rmse_wavg)
})

test_that("cprc_class adapts to abrupt hyperplane drift and returns valid probabilities", {
  set.seed(42)
  dat <- make_seq_plane(n1 = 120, n2 = 200, blocks = 8, flippoint = c(5),
                        dimen = 6, delta = c(0.02, 2.0), flipprob = 0.05)
  mods <- lapply(dat$train_batches[1:6], function(df) {
    suppressWarnings(glm(y ~ ., data = df, family = binomial()))
  })
  Fm_tr <- sapply(mods, function(m) predict(m, newdata = dat$train_batches[[7]], type = "link"))
  Fm_te <- sapply(mods, function(m) predict(m, newdata = dat$test_batches[[8]], type = "link"))
  y_tr <- dat$train_batches[[7]]$y

  fit <- cprc_class(y = y_tr, Fm = Fm_tr, lambda = 5, penalty = "lasso", cv = FALSE)
  expect_s3_class(fit, "cprc_class")
  expect_equal(sum(fit$w), 1, tolerance = 1e-6)
  expect_true(all(fit$w >= 0))
  expect_true(mean(fit$w[5:6]) > mean(fit$w[1:4]))

  probs <- predict(fit, newdata = Fm_te, type = "prob")
  cls   <- predict(fit, newdata = Fm_te, type = "class")
  expect_true(all(probs >= 0 & probs <= 1))
  expect_true(all(cls %in% c(0, 1)))
})

test_that("cprc_unsupervised converges, down-weights outliers, and detects active voxels", {
  spm_path <- system.file("extdata", "synth_spm.csv", package = "cprc")
  skip_if(spm_path == "", "synth_spm.csv not found")
  X <- as.matrix(read.csv(spm_path))

  fit <- cprc_unsupervised(X, lambda = 20, max_iter = 300)
  expect_s3_class(fit, "cprc_unsupervised")
  expect_true(fit$converged)
  expect_equal(sum(fit$w), 1, tolerance = 1e-6)
  expect_true(mean(fit$w[18:20]) < 0.025)
  expect_true(mean(fit$w[1:17]) > 0.052)

  set.seed(123)
  bt <- cprc_boot_test(X, w = fit, B = 100, alpha = 0.05, adjust_method = "holm")
  expect_equal(length(bt$active_voxels), nrow(X))
  expect_true(sum(bt$active_voxels[1:43]) >= 35L)
  expect_true(sum(bt$active_voxels[44:121]) <= 3L)
})

test_that("make_seq_sea produces valid binary streaming batches", {
  set.seed(7)
  sea <- make_seq_sea(batch_size = 40, num_batches = 4)
  expect_equal(length(sea$train_batches), 4L)
  expect_equal(nrow(sea$train_batches[[1]]), 40L)
  expect_true(all(sea$train_batches[[1]]$y %in% c(0, 1)))
})
