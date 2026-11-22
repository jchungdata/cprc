#' Constrained Penalized Regression Combiner for Streaming Regression
#'
#' Fits a model combiner for streaming regression under concept drift with non-negativity
#' and sum-to-one constraints on the weights.
#'
#' @param y Numeric response vector of length \eqn{n}.
#' @param Fm Numeric matrix of dimension \eqn{n \times M} containing predictions from \eqn{M} base models.
#' @param lambda Regularization parameter \eqn{\lambda \ge 0}. If \code{NULL}, tuned via K-fold cross-validation.
#' @param penalty Character specifying combiner penalty: \code{"lasso"} (MC.Lasso1+, default),
#'   \code{"ridge"} (MC.Ridge1+), \code{"lse"} (MC.Lse1+), or \code{"wavg"} (MC.WAvg).
#' @param type Character specifying penalty structure: \code{"pairwise"} (\eqn{\|\mathbf{D}\mathbf{w}\|}) or
#'   \code{"identity"} (\eqn{\|\mathbf{w}\|}). Default is \code{"pairwise"}.
#' @param cv Logical. Whether to perform cross-validation to select \eqn{\lambda}. Default is \code{TRUE} if \code{lambda} is NULL.
#' @param lambda_seq Optional grid of \eqn{\lambda} candidate values for cross-validation.
#' @param folds Number of folds for cross-validation. Default is 10.
#' @param cv_rule Selection rule for cross-validation: \code{"1se"} (default) or \code{"min"}.
#' @param rho Augmented Lagrangian parameter for ADMM. Default is 100.
#' @param max_step Maximum iterations for ADMM solver. Default is 5000.
#' @param conv_eps Convergence threshold. Default is 1e-8.
#' @param verbose Logical. Whether to print progress messages. Default is \code{FALSE}.
#'
#' @return An object of class \code{"cprc_reg"} containing:
#' \item{w}{Estimated combiner weight vector of length \eqn{M}.}
#' \item{lambda}{Selected regularization parameter \eqn{\lambda}.}
#' \item{penalty}{Penalty type utilized.}
#' \item{type}{Penalty structure (\code{"pairwise"} or \code{"identity"}).}
#' \item{angle}{Angular departure from uniform weighting \eqn{\eta(\mathbf{w}^*, \hat{\mathbf{w}})} in degrees.}
#' \item{clusters}{Identified batch clusters and membership indices.}
#' \item{fitted_values}{Fitted values on the training data \eqn{\hat{\mathbf{y}} = \mathbf{F}\hat{\mathbf{w}}}.}
#' \item{residuals}{Residual vector \eqn{\mathbf{y} - \hat{\mathbf{y}}}.}
#' \item{call}{Matched function call.}
#'
#' @references
#' Wang, L.-Y., Park, C., Yeon, K., and Choi, H. (2017). Tracking concept drift using a constrained
#' penalized regression combiner. \emph{Computational Statistics & Data Analysis}, 108, 52--69.
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' dat <- make_seq_reg(n = 100, M = 11, change_points = c(6))
#' mods <- lapply(dat$seqdata[1:10], function(df) lm(y ~ ., data = df))
#' Fm <- sapply(mods, function(m) predict(m, newdata = dat$seqdata[[11]]))
#' fit <- cprc_reg(y = dat$seqdata[[11]]$y, Fm = Fm, penalty = "lasso")
#' print(fit)
#' plot(fit)
#' }
cprc_reg <- function(y, Fm, lambda = NULL, penalty = c("lasso", "ridge", "lse", "wavg"),
                     type = c("pairwise", "identity"), cv = is.null(lambda),
                     lambda_seq = NULL, folds = 10, cv_rule = c("1se", "min"), rho = 100,
                     max_step = 5000, conv_eps = 1e-8, verbose = FALSE) {
  cl <- match.call()
  penalty <- match.arg(penalty)
  type <- match.arg(type)
  cv_rule <- match.arg(cv_rule)

  Fm <- as.matrix(Fm)
  y <- as.vector(y)
  n <- nrow(Fm)
  M <- ncol(Fm)

  if (length(y) != n) {
    stop("Length of y must match number of rows in Fm.")
  }

  cv_res <- NULL

  # 1. Unpenalized Constrained Least Squares (MC.Lse1+)
  if (penalty == "lse") {
    sol <- lse_qp(y, Fm)
    w <- sol$w
    lambda_val <- 0
    fit_clusters <- list(membership = seq_len(M), csize = rep(1, M), no = M)

  # 2. Accuracy-based Weighted Average (MC.WAvg)
  } else if (penalty == "wavg") {
    sol <- weighted_avg(y, Fm)
    w <- sol$w
    lambda_val <- 0
    fit_clusters <- list(membership = seq_len(M), csize = rep(1, M), no = M)

  # 3. Constrained Ridge Combiner (MC.Ridge1+)
  } else if (penalty == "ridge") {
    if (cv && is.null(lambda)) {
      cv_res <- cv_cprc_reg(y, Fm, lambda_seq = lambda_seq, folds = folds,
                            penalty = "ridge", type = type, cv_rule = cv_rule, rho = rho,
                            conv_eps = conv_eps, verbose = verbose)
      lambda_val <- cv_res$best_lambda
    } else {
      lambda_val <- if (is.null(lambda)) 1 else lambda
    }

    if (type == "identity") {
      sol <- ridge_qp(y, Fm, lambda = lambda_val)
      w <- sol$w
      fit_clusters <- list(membership = seq_len(M), csize = rep(1, M), no = M)
    } else {
      sol <- admm_cprc(y, Fm, lambda = lambda_val, method = "regression",
                       norm = "l2", rho = rho, max_step = max_step,
                       conv_eps = conv_eps, verbose = verbose)
      w <- sol$w
      fit_clusters <- sol$clusters
    }

  # 4. Generalized Fused Lasso Combiner (MC.Lasso1+)
  } else if (penalty == "lasso") {
    if (cv && is.null(lambda)) {
      cv_res <- cv_cprc_reg(y, Fm, lambda_seq = lambda_seq, folds = folds,
                            penalty = "lasso", type = type, cv_rule = cv_rule, rho = rho,
                            conv_eps = conv_eps, verbose = verbose)
      lambda_val <- cv_res$best_lambda
    } else {
      lambda_val <- if (is.null(lambda)) 1 else lambda
    }

    sol <- admm_cprc(y, Fm, lambda = lambda_val, method = "regression",
                     norm = "l1", rho = rho, max_step = max_step,
                     conv_eps = conv_eps, verbose = verbose)
    w <- sol$w
    fit_clusters <- sol$clusters
  }

  drift_ang <- concept_drift_angle(w, M = M, unit = "degrees")
  fitted_val <- drop(Fm %*% w)
  resid_val <- y - fitted_val

  out <- list(
    w = w,
    lambda = lambda_val,
    penalty = penalty,
    type = type,
    angle = drift_ang,
    clusters = fit_clusters,
    fitted_values = fitted_val,
    residuals = resid_val,
    cv = cv_res,
    n = n,
    M = M,
    call = cl
  )
  class(out) <- "cprc_reg"
  return(out)
}
