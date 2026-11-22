#' Bootstrap-Based Weighted t-Test for Statistical Parametric Maps
#'
#' Performs voxel-wise significance testing on aggregated statistical parametric maps using subject
#' weights estimated by \code{\link{cprc_unsupervised}}. Computes the weighted \eqn{t}-statistic and
#' evaluates statistical significance via a parametric bootstrap under the empirical null distribution.
#'
#' @param X Numeric matrix of dimension \eqn{M \times N}, containing regression coefficients or statistics
#'   at \eqn{M} voxels across \eqn{N} subjects.
#' @param w Numeric vector of length \eqn{N} containing subject aggregation weights. If an object of class
#'   \code{"cprc_unsupervised"} is provided, weights are automatically extracted.
#' @param B Integer. Number of bootstrap replications per voxel. Default is 100.
#' @param alpha Significance level for identifying active voxels. Default is 0.05.
#' @param adjust_method Character string specifying the multiple testing adjustment method:
#'   \code{"holm"} (default), \code{"fdr"} (Benjamini-Hochberg), or \code{"none"}.
#'
#' @return A list containing:
#' \item{t_stat}{Numeric vector of weighted \eqn{t}-statistics for each voxel.}
#' \item{p_values}{Raw empirical bootstrap \eqn{p}-values for each voxel.}
#' \item{p_adjusted}{Adjusted \eqn{p}-values according to \code{adjust_method}.}
#' \item{active_voxels}{Logical vector indicating statistically significant voxels at level \code{alpha}.}
#' \item{weighted_mean}{Weighted mean coefficient map \eqn{\hat{\boldsymbol{\beta}}^w}.}
#' \item{weighted_sd}{Weighted standard deviation map \eqn{\mathbf{S}^w}.}
#'
#' @references
#' Wang, L.-Y., Chung, J., Park, C., Choi, H., Rodrigue, A. L., Pierce, J. E., Clementz, B. A.,
#' and McDowell, J. E. (2018). Regularized aggregation of statistical parametric maps.
#' \emph{Human Brain Mapping}, 39(12), 4811--4825.
#'
#' @export
cprc_boot_test <- function(X, w, B = 100, alpha = 0.05,
                           adjust_method = c("holm", "fdr", "none")) {
  adjust_method <- match.arg(adjust_method)

  if (inherits(w, "cprc_unsupervised")) {
    w <- w$w
  }
  X <- as.matrix(X)
  M <- nrow(X)
  N <- ncol(X)

  if (length(w) != N) {
    stop("Length of weight vector w must match the number of columns in X.")
  }

  w <- w / sum(w)

  # Compute weighted mean and weighted standard deviation for each voxel
  weighted_mean <- drop(X %*% w)
  weighted_sd <- numeric(M)

  for (v in seq_len(M)) {
    weighted_sd[v] <- sqrt(sum(w * (X[v, ] - weighted_mean[v])^2))
  }
  weighted_sd[weighted_sd <= 0] <- 1e-6

  # Compute observed weighted t-statistics
  t_stat <- weighted_mean / (weighted_sd / sqrt(N))

  # Parametric bootstrap under the null H0: beta_v = 0
  p_values <- numeric(M)

  for (v in seq_len(M)) {
    obs_t <- t_stat[v]
    sigma_v <- weighted_sd[v]

    boot_t <- numeric(B)
    for (b in seq_len(B)) {
      b_sample <- stats::rnorm(N, mean = 0, sd = sigma_v)
      boot_t[b] <- stats::t.test(b_sample)$statistic
    }

    prop_greater <- mean(boot_t >= obs_t)
    prop_lesser <- mean(boot_t <= obs_t)
    p_values[v] <- min(1, 2 * min(prop_greater, prop_lesser))
  }

  if (adjust_method == "none") {
    p_adjusted <- p_values
  } else {
    p_adjusted <- stats::p.adjust(p_values, method = adjust_method)
  }

  active_voxels <- (p_adjusted < alpha)

  return(list(
    t_stat = t_stat,
    p_values = p_values,
    p_adjusted = p_adjusted,
    active_voxels = active_voxels,
    weighted_mean = weighted_mean,
    weighted_sd = weighted_sd,
    adjust_method = adjust_method,
    alpha = alpha
  ))
}
