#' Constrained Penalized Regression Combiner for Streaming Classification
#'
#' Fits a model combiner for binary classification in streaming environments with concept drift,
#' optimizing the negative Bernoulli log-likelihood loss under simplex constraints on the weights.
#'
#' @param y Binary response vector with values in \{0, 1\} of length \eqn{n}.
#' @param Fm Numeric matrix of dimension \eqn{n \times M} containing predicted probabilities
#'   or continuous decision scores from \eqn{M} base classifiers.
#' @param lambda Regularization parameter \eqn{\lambda \ge 0}. If \code{NULL}, tuned via K-fold cross-validation.
#' @param penalty Character specifying combiner penalty: \code{"lasso"} (fused lasso difference penalty, default)
#'   or \code{"ridge"} (pairwise \eqn{\ell_2} difference penalty).
#' @param cv Logical. Whether to perform cross-validation to select \eqn{\lambda}. Default is \code{TRUE} if \code{lambda} is NULL.
#' @param lambda_seq Optional grid of \eqn{\lambda} candidate values for cross-validation.
#' @param folds Number of folds for cross-validation. Default is 10.
#' @param cv_rule Selection rule for cross-validation: \code{"1se"} (default) or \code{"min"}.
#' @param rho Augmented Lagrangian parameter for ADMM. Default is 100.
#' @param max_step Maximum iterations for ADMM solver. Default is 5000.
#' @param conv_eps Convergence threshold. Default is 1e-8.
#' @param verbose Logical. Whether to print progress messages. Default is \code{FALSE}.
#'
#' @return An object of class \code{"cprc_class"} containing:
#' \item{w}{Estimated combiner weight vector of length \eqn{M}.}
#' \item{lambda}{Selected regularization parameter \eqn{\lambda}.}
#' \item{penalty}{Penalty type utilized.}
#' \item{angle}{Angular departure from uniform weighting \eqn{\eta(\mathbf{w}^*, \hat{\mathbf{w}})} in degrees.}
#' \item{clusters}{Identified batch clusters and membership indices.}
#' \item{fitted_probabilities}{Predicted event probabilities \eqn{\hat{p}_i = 1 / (1 + \exp(-\mathbf{F}_i \hat{\mathbf{w}}))}.}
#' \item{fitted_classes}{Predicted class labels (0 or 1) based on a 0.5 probability cutoff.}
#' \item{call}{Matched function call.}
#'
#' @references
#' Wang, L.-Y., Park, C., Choi, H., and Yeon, K. (2016). A classifier ensemble for concept drift
#' using a constrained penalized regression combiner. \emph{Procedia Computer Science}, 91, 252--259.
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(42)
#' dat <- make_seq_plane(n1 = 100, n2 = 200, blocks = 11, flippoint = c(6))
#' mods <- lapply(dat$train_batches[1:10], function(df) glm(y ~ ., data = df, family = binomial()))
#' Fm <- sapply(mods, function(m) predict(m, newdata = dat$train_batches[[11]], type = "link"))
#' fit <- cprc_class(y = dat$train_batches[[11]]$y, Fm = Fm, penalty = "lasso")
#' print(fit)
#' }
cprc_class <- function(y, Fm, lambda = NULL, penalty = c("lasso", "ridge"),
                       cv = is.null(lambda), lambda_seq = NULL, folds = 10, cv_rule = c("1se", "min"),
                       rho = 100, max_step = 5000, conv_eps = 1e-8, verbose = FALSE) {
  cl <- match.call()
  penalty <- match.arg(penalty)
  cv_rule <- match.arg(cv_rule)

  Fm <- as.matrix(Fm)
  y <- as.vector(y)
  n <- nrow(Fm)
  M <- ncol(Fm)

  if (length(y) != n) {
    stop("Length of y must match number of rows in Fm.")
  }
  if (!all(y %in% c(0, 1))) {
    stop("Classification response y must be binary with values in {0, 1}.")
  }

  norm_type <- if (penalty == "lasso") "l1" else "l2"
  cv_res <- NULL

  if (cv && is.null(lambda)) {
    cv_res <- cv_cprc_class(y, Fm, lambda_seq = lambda_seq, folds = folds,
                            penalty = penalty, cv_rule = cv_rule, rho = rho, conv_eps = conv_eps,
                            verbose = verbose)
    lambda_val <- cv_res$best_lambda
  } else {
    lambda_val <- if (is.null(lambda)) 1 else lambda
  }

  sol <- admm_cprc(y, Fm, lambda = lambda_val, method = "classification",
                   norm = norm_type, rho = rho, max_step = max_step,
                   conv_eps = conv_eps, verbose = verbose)
  w <- sol$w

  h <- drop(Fm %*% w)
  prob <- 1 / (1 + exp(-pmax(pmin(h, 50), -50)))
  pred_class <- ifelse(prob >= 0.5, 1, 0)
  drift_ang <- concept_drift_angle(w, M = M, unit = "degrees")

  out <- list(
    w = w,
    lambda = lambda_val,
    penalty = penalty,
    angle = drift_ang,
    clusters = sol$clusters,
    fitted_probabilities = prob,
    fitted_classes = pred_class,
    cv = cv_res,
    n = n,
    M = M,
    call = cl
  )
  class(out) <- "cprc_class"
  return(out)
}
