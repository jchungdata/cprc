#' K-Fold Cross-Validation Routines for cprc
#'
#' Selects the optimal regularization parameter \eqn{\lambda} via K-fold cross-validation
#' for supervised streaming regression, supervised streaming classification, and
#' quadratic programming ridge combiners.
#'
#' @name cv_cprc
NULL

#' Cross-Validation for Streaming Regression Combiners
#'
#' @param y Numeric response vector of length \eqn{n}.
#' @param Fm Numeric matrix of dimension \eqn{n \times M} of base learner predictions.
#' @param lambda_seq Optional numeric vector of candidate \eqn{\lambda} values.
#' @param folds Integer number of CV folds. Default is 10.
#' @param penalty Character: \code{"lasso"} or \code{"ridge"}.
#' @param type Character: \code{"pairwise"} or \code{"identity"}.
#' @param cv_rule Selection rule for \eqn{\lambda}: \code{"1se"} (largest \eqn{\lambda} within
#'   one standard error of the minimum CV error, default) or \code{"min"} (minimum CV error).
#' @param rho Augmented Lagrangian parameter for ADMM. Default is 100.
#' @param conv_eps Convergence tolerance. Default is 1e-8.
#' @param verbose Logical. Print progress if \code{TRUE}.
#' @return A list with \code{best_lambda}, \code{lambda_min}, \code{lambda_1se},
#'   \code{lambda_seq}, \code{cv_error}, \code{cv_se}, and \code{cv_rule}.
#' @export
cv_cprc_reg <- function(y, Fm, lambda_seq = NULL, folds = 10,
                        penalty = c("lasso", "ridge"),
                        type = c("pairwise", "identity"),
                        cv_rule = c("1se", "min"),
                        rho = 100, conv_eps = 1e-8, verbose = FALSE) {
  penalty <- match.arg(penalty)
  type <- match.arg(type)
  cv_rule <- match.arg(cv_rule)
  Fm <- as.matrix(Fm)
  y <- as.vector(y)
  n <- nrow(Fm)

  if (is.null(lambda_seq)) {
    lambda_seq <- c(0.1, 0.5, 1, 5, 10, 25, 50, 100, 250, 500)
  } else {
    lambda_seq <- sort(unique(as.numeric(lambda_seq)))
  }

  fold_id <- cv_folds(n, K = folds)
  K <- max(fold_id)
  cv_errors <- numeric(length(lambda_seq))
  cv_se <- numeric(length(lambda_seq))

  for (j in seq_along(lambda_seq)) {
    lam <- lambda_seq[j]
    fold_err <- numeric(K)

    for (k in seq_len(K)) {
      val_idx <- which(fold_id == k)
      tr_idx <- which(fold_id != k)

      if (penalty == "ridge" && type == "identity") {
        sol <- ridge_qp(y[tr_idx], Fm[tr_idx, , drop = FALSE], lambda = lam)
      } else {
        norm_type <- if (penalty == "lasso") "l1" else "l2"
        sol <- admm_cprc(y[tr_idx], Fm[tr_idx, , drop = FALSE], lambda = lam,
                         method = "regression", norm = norm_type,
                         rho = rho, conv_eps = conv_eps, verbose = FALSE)
      }
      pred <- drop(Fm[val_idx, , drop = FALSE] %*% sol$w)
      fold_err[k] <- calc_rmse(pred, y[val_idx])
    }
    cv_errors[j] <- mean(fold_err)
    cv_se[j] <- if (K > 1L) stats::sd(fold_err) / sqrt(K) else 0
    if (verbose) message(sprintf("lambda = %.4f | CV RMSE = %.4f (SE = %.4f)", lam, cv_errors[j], cv_se[j]))
  }

  best_idx_min <- which.min(cv_errors)
  lambda_min <- lambda_seq[best_idx_min]
  se_thresh <- cv_errors[best_idx_min] + cv_se[best_idx_min]
  best_idx_1se <- max(which(cv_errors <= se_thresh))
  lambda_1se <- lambda_seq[best_idx_1se]
  best_lambda <- if (cv_rule == "1se") lambda_1se else lambda_min

  return(list(
    best_lambda = best_lambda,
    lambda_min = lambda_min,
    lambda_1se = lambda_1se,
    lambda_seq = lambda_seq,
    cv_error = cv_errors,
    cv_se = cv_se,
    cv_rule = cv_rule
  ))
}

#' Cross-Validation for Streaming Classification Combiners
#'
#' @param y Binary response vector in \{0, 1\} of length \eqn{n}.
#' @param Fm Numeric matrix of dimension \eqn{n \times M} of base classifier scores.
#' @param lambda_seq Optional numeric vector of candidate \eqn{\lambda} values.
#' @param folds Integer number of CV folds. Default is 10.
#' @param penalty Character: \code{"lasso"} or \code{"ridge"}.
#' @param cv_rule Selection rule for \eqn{\lambda}: \code{"1se"} (largest \eqn{\lambda} within
#'   one standard error of the minimum CV error, default) or \code{"min"} (minimum CV error).
#' @param rho Augmented Lagrangian parameter for ADMM. Default is 100.
#' @param conv_eps Convergence tolerance. Default is 1e-8.
#' @param verbose Logical. Print progress if \code{TRUE}.
#' @return A list with \code{best_lambda}, \code{lambda_min}, \code{lambda_1se},
#'   \code{lambda_seq}, \code{cv_error}, \code{cv_se}, and \code{cv_rule}.
#' @export
cv_cprc_class <- function(y, Fm, lambda_seq = NULL, folds = 10,
                          penalty = c("lasso", "ridge"),
                          cv_rule = c("1se", "min"),
                          rho = 100, conv_eps = 1e-8, verbose = FALSE) {
  penalty <- match.arg(penalty)
  cv_rule <- match.arg(cv_rule)
  Fm <- as.matrix(Fm)
  y <- as.vector(y)
  n <- nrow(Fm)

  if (is.null(lambda_seq)) {
    lambda_seq <- c(0.1, 0.5, 1, 2, 5, 10, 25, 50, 100)
  } else {
    lambda_seq <- sort(unique(as.numeric(lambda_seq)))
  }

  norm_type <- if (penalty == "lasso") "l1" else "l2"
  fold_id <- cv_folds(n, K = folds)
  K <- max(fold_id)
  cv_errors <- numeric(length(lambda_seq))
  cv_se <- numeric(length(lambda_seq))

  for (j in seq_along(lambda_seq)) {
    lam <- lambda_seq[j]
    fold_err <- numeric(K)

    for (k in seq_len(K)) {
      val_idx <- which(fold_id == k)
      tr_idx <- which(fold_id != k)

      sol <- admm_cprc(y[tr_idx], Fm[tr_idx, , drop = FALSE], lambda = lam,
                       method = "classification", norm = norm_type,
                       rho = rho, conv_eps = conv_eps, verbose = FALSE)
      h <- drop(Fm[val_idx, , drop = FALSE] %*% sol$w)
      prob <- 1 / (1 + exp(-pmax(pmin(h, 50), -50)))
      prob_clip <- pmin(pmax(prob, 1e-12), 1 - 1e-12)
      fold_err[k] <- -mean(y[val_idx] * log(prob_clip) + (1 - y[val_idx]) * log(1 - prob_clip))
    }
    cv_errors[j] <- mean(fold_err)
    cv_se[j] <- if (K > 1L) stats::sd(fold_err) / sqrt(K) else 0
    if (verbose) message(sprintf("lambda = %.4f | CV LogLoss = %.4f (SE = %.4f)", lam, cv_errors[j], cv_se[j]))
  }

  best_idx_min <- which.min(cv_errors)
  lambda_min <- lambda_seq[best_idx_min]
  se_thresh <- cv_errors[best_idx_min] + cv_se[best_idx_min]
  best_idx_1se <- max(which(cv_errors <= se_thresh))
  lambda_1se <- lambda_seq[best_idx_1se]
  best_lambda <- if (cv_rule == "1se") lambda_1se else lambda_min

  return(list(
    best_lambda = best_lambda,
    lambda_min = lambda_min,
    lambda_1se = lambda_1se,
    lambda_seq = lambda_seq,
    cv_error = cv_errors,
    cv_se = cv_se,
    cv_rule = cv_rule
  ))
}

#' Cross-Validation for Quadratic Programming Ridge Combiner
#'
#' @param y Numeric response vector of length \eqn{n}.
#' @param Fm Numeric matrix of dimension \eqn{n \times M}.
#' @param folds Integer number of CV folds. Default is 10.
#' @param lambda_seq Numeric vector of candidate \eqn{\lambda} values.
#' @param cv_rule Selection rule for \eqn{\lambda}: \code{"1se"} (default) or \code{"min"}.
#' @return A list with \code{best_lambda}, \code{lambda_min}, \code{lambda_1se},
#'   \code{lambda_seq}, \code{cv_error}, \code{cv_se}, and \code{cv_rule}.
#' @export
cv_cprc_qp <- function(y, Fm, folds = 10, lambda_seq = 2 * 10^seq(-3, 3, by = 0.5),
                       cv_rule = c("1se", "min")) {
  cv_rule <- match.arg(cv_rule)
  Fm <- as.matrix(Fm)
  y <- as.vector(y)
  n <- nrow(Fm)
  lambda_seq <- sort(unique(as.numeric(lambda_seq)))
  fold_id <- cv_folds(n, K = folds)
  K <- max(fold_id)
  cv_errors <- numeric(length(lambda_seq))
  cv_se <- numeric(length(lambda_seq))

  for (j in seq_along(lambda_seq)) {
    lam <- lambda_seq[j]
    fold_err <- numeric(K)
    for (k in seq_len(K)) {
      val_idx <- which(fold_id == k)
      tr_idx <- which(fold_id != k)
      sol <- ridge_qp(y[tr_idx], Fm[tr_idx, , drop = FALSE], lambda = lam)
      pred <- drop(Fm[val_idx, , drop = FALSE] %*% sol$w)
      fold_err[k] <- mean((y[val_idx] - pred)^2)
    }
    cv_errors[j] <- mean(fold_err)
    cv_se[j] <- if (K > 1L) stats::sd(fold_err) / sqrt(K) else 0
  }

  best_idx_min <- which.min(cv_errors)
  lambda_min <- lambda_seq[best_idx_min]
  se_thresh <- cv_errors[best_idx_min] + cv_se[best_idx_min]
  best_idx_1se <- max(which(cv_errors <= se_thresh))
  lambda_1se <- lambda_seq[best_idx_1se]
  best_lambda <- if (cv_rule == "1se") lambda_1se else lambda_min

  return(list(
    best_lambda = best_lambda,
    lambda_min = lambda_min,
    lambda_1se = lambda_1se,
    lambda_seq = lambda_seq,
    cv_error = cv_errors,
    cv_se = cv_se,
    cv_rule = cv_rule
  ))
}
