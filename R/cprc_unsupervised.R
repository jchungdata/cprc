#' Regularized Unsupervised Aggregation of Statistical Parametric Maps
#'
#' Aggregates high-dimensional parameter maps (such as voxel-wise SPMs in neuroimaging) across
#' multiple subjects when the ground-truth target map is unobserved. Employs constrained ridge
#' regularization on the probability simplex coupled with iterative alternating optimal scoring
#' initialized by the median map.
#'
#' @param X Numeric matrix of dimension \eqn{M \times N}, where \eqn{M} is the number of features/voxels
#'   and \eqn{N} is the number of subjects or observation sources.
#' @param lambda Regularization parameter \eqn{\lambda \ge 0}. If \code{NULL}, tuned via K-fold cross-validation.
#' @param lambda_seq Grid of candidate \eqn{\lambda} values. Default is \code{2 * 10^seq(-3, 3, by = 0.5)}.
#' @param max_iter Maximum alternating iterations. Default is 1000.
#' @param tol Relative convergence tolerance. Default is 1e-6.
#' @param warmup Minimum number of iterations before assessing convergence. Default is 20.
#' @param folds Number of folds for cross-validation when tuning \eqn{\lambda}. Default is 10.
#' @param verbose Logical. Whether to print iteration progress. Default is \code{FALSE}.
#'
#' @return An object of class \code{"cprc_unsupervised"} containing:
#' \item{w}{Estimated optimal subject weight vector of length \eqn{N}.}
#' \item{combined_map}{Aggregated representative map \eqn{\hat{\mathbf{y}} = \mathbf{X}\hat{\mathbf{w}}} of length \eqn{M}.}
#' \item{w_history}{Matrix recording estimated weights across alternating iterations.}
#' \item{obj_history}{Vector of objective function values across iterations.}
#' \item{lambda}{Selected regularization parameter \eqn{\lambda}.}
#' \item{converged}{Logical indicating whether the algorithm converged within \code{max_iter}.}
#' \item{iterations}{Total number of alternating iterations completed.}
#' \item{call}{Matched function call.}
#'
#' @references
#' Wang, L.-Y., Chung, J., Park, C., Choi, H., Rodrigue, A. L., Pierce, J. E., Clementz, B. A.,
#' and McDowell, J. E. (2018). Regularized aggregation of statistical parametric maps.
#' \emph{Human Brain Mapping}, 39(12), 4811--4825.
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(121 * 20), nrow = 121, ncol = 20)
#' fit_un <- cprc_unsupervised(X, lambda = 20, max_iter = 50)
#' print(fit_un)
#' plot(fit_un, type = "elbow")
#' }
cprc_unsupervised <- function(X, lambda = NULL, lambda_seq = 2 * 10^seq(-3, 3, by = 0.5),
                              max_iter = 1000, tol = 1e-6, warmup = 20,
                              folds = 10, verbose = FALSE) {
  cl <- match.call()
  X <- as.matrix(X)
  M <- nrow(X)
  N <- ncol(X)

  if (N < 2L) {
    stop("X must have at least 2 subjects (columns) to perform aggregation.")
  }

  # Step 1: Initialize y with the voxel-wise median vector of X and fix its L2 scale
  # for optimal scoring so the latent map does not degenerate to uniform weighting
  y_prev <- apply(X, 1, stats::median)
  y_scale <- sqrt(sum(y_prev^2))
  w <- rep(1 / N, N)

  # Tuning lambda
  tune_lambda <- is.null(lambda)
  if (tune_lambda) {
    cv_res <- cv_cprc_qp(y = y_prev, Fm = X, folds = folds, lambda_seq = lambda_seq)
    curr_lambda <- cv_res$best_lambda
  } else {
    curr_lambda <- lambda
  }

  fit_init <- ridge_qp(y = y_prev, Fm = X, lambda = curr_lambda)
  w <- fit_init$w
  xw <- drop(X %*% w)
  xw_norm <- sqrt(sum(xw^2))
  y <- if (y_scale > 0 && xw_norm > 0) xw * (y_scale / xw_norm) else xw

  w_store <- matrix(NA_real_, nrow = N, ncol = max_iter + 1L)
  w_store[, 1L] <- w
  obj_history <- numeric(max_iter + 1L)
  obj_history[1L] <- sum((y - xw)^2) + curr_lambda * sum(w^2)

  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(max_iter)) {
    if (verbose && (iter %% 10L == 0L || iter <= 5L)) {
      message(sprintf("cprc_unsupervised iteration: %d", iter))
    }

    # Step 2: Update lambda during initial iterations if requested
    if (tune_lambda && iter <= 5L) {
      cv_res <- cv_cprc_qp(y = y, Fm = X, folds = folds, lambda_seq = lambda_seq)
      curr_lambda <- cv_res$best_lambda
    }

    # Step 3: Solve for weight w given current latent map y
    fit_step <- ridge_qp(y = y, Fm = X, lambda = curr_lambda)
    w_new <- fit_step$w

    # Step 4: Update norm-constrained latent map y^{(t)} and evaluate optimal scoring objective
    xw_new <- drop(X %*% w_new)
    xw_new_norm <- sqrt(sum(xw_new^2))
    y_new <- if (y_scale > 0 && xw_new_norm > 0) xw_new * (y_scale / xw_new_norm) else xw_new
    curr_obj <- sum((y_new - xw_new)^2) + curr_lambda * sum(w_new^2)

    w_store[, iter + 1L] <- w_new
    obj_history[iter + 1L] <- curr_obj

    # Check stopping criterion after warmup
    if (iter >= warmup) {
      prev_obj <- obj_history[iter]
      rel_diff <- abs(curr_obj - prev_obj) / (prev_obj + 1e-12)
      w_diff <- max(abs(w_new - w))
      if (rel_diff < tol && w_diff < tol) {
        converged <- TRUE
        w_store <- w_store[, seq_len(iter + 1L), drop = FALSE]
        obj_history <- obj_history[seq_len(iter + 1L)]
        if (verbose) {
          message(sprintf("Convergence reached at iteration %d (rel_diff = %.2e)", iter, rel_diff))
        }
        break
      }
    }

    w <- w_new
    y <- y_new
  }

  if (!converged) {
    w_store <- w_store[, seq_len(iter + 1L), drop = FALSE]
    obj_history <- obj_history[seq_len(iter + 1L)]
  }

  final_w <- w_store[, ncol(w_store)]
  final_map <- drop(X %*% final_w)

  structure(
    list(
      w = final_w,
      combined_map = final_map,
      w_history = w_store,
      obj_history = obj_history,
      lambda = curr_lambda,
      converged = converged,
      iterations = iter,
      M = M,
      N = N,
      call = cl
    ),
    class = "cprc_unsupervised"
  )
}
