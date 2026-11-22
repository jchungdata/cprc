#' Quadratic Programming Solvers for Ridge and LSE Combiners
#'
#' Implements constrained ridge regression combiners (identity and pairwise difference penalties),
#' constrained least squares, and accuracy-based weighted average combiners on the probability simplex.
#'
#' @name qp_solvers
NULL

#' Constrained Ridge Combiner via Quadratic Programming (Identity Penalty)
#'
#' Solves the constrained \eqn{\ell_2} penalized regression combiner on the unit simplex:
#' \deqn{\min_{\mathbf{w}} \|\mathbf{y} - \mathbf{F}\mathbf{w}\|_2^2 + \lambda \|\mathbf{w}\|_2^2 \quad \text{s.t.} \quad \sum_{m=1}^M w_m = 1, \; w_m \ge 0}
#' using \code{\link[quadprog]{solve.QP}}.
#'
#' @param y Numeric vector of responses of length \eqn{n}.
#' @param Fm Numeric matrix of base learner predictions of dimension \eqn{n \times M}.
#' @param lambda Regularization tuning parameter \eqn{\lambda \ge 0}.
#' @return A list containing \code{w} (estimated weight vector).
#' @references
#' Yeon, K., Song, M. S., Kim, Y., Choi, H., and Park, C. (2010). Model averaging via penalized regression
#' for tracking concept drift. \emph{Journal of Computational and Graphical Statistics}, 19(2), 457--473.
#' @export
ridge_qp <- function(y, Fm, lambda) {
  Fm <- as.matrix(Fm)
  y <- as.vector(y)
  m <- ncol(Fm)

  dvec <- drop(crossprod(Fm, y))
  bvec <- c(1, rep(0, m))
  Amat <- rbind(rep(1, m), diag(m))
  Dmat <- crossprod(Fm) + lambda * diag(m)

  # Ensure Dmat is strictly positive definite for solve.QP
  mineig <- min(eigen(Dmat, symmetric = TRUE, only.values = TRUE)$values)
  if (mineig <= 0) {
    Dmat <- Dmat + (abs(mineig) + 1e-6) * diag(m)
  }

  res <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq = 1, factorized = FALSE)
  w <- res$solution
  w[w < 0] <- 0
  w <- w / sum(w)

  return(list(w = w))
}

#' Constrained Ridge Combiner with Pairwise Difference Penalty
#'
#' Solves the constrained pairwise difference \eqn{\ell_2} combiner on the unit simplex:
#' \deqn{\min_{\mathbf{w}} \|\mathbf{y} - \mathbf{F}\mathbf{w}\|_2^2 + \lambda \|\mathbf{D}\mathbf{w}\|_2^2 \quad \text{s.t.} \quad \sum_{m=1}^M w_m = 1, \; w_m \ge 0}
#' using \code{\link[quadprog]{solve.QP}}.
#'
#' @param y Numeric vector of responses of length \eqn{n}.
#' @param Fm Numeric matrix of base learner predictions of dimension \eqn{n \times M}.
#' @param lambda Regularization tuning parameter \eqn{\lambda \ge 0}.
#' @param D Difference matrix. If \code{NULL}, generated via \code{\link{diff_mat}}.
#' @return A list containing \code{w} (estimated weight vector).
#' @export
ridge_qp_diff <- function(y, Fm, lambda, D = NULL) {
  Fm <- as.matrix(Fm)
  y <- as.vector(y)
  m <- ncol(Fm)

  if (is.null(D)) {
    D <- diff_mat(m)
  }

  dvec <- drop(crossprod(Fm, y))
  bvec <- c(1, rep(0, m))
  Amat <- rbind(rep(1, m), diag(m))
  Dmat <- crossprod(Fm) + lambda * crossprod(D)

  mineig <- min(eigen(Dmat, symmetric = TRUE, only.values = TRUE)$values)
  if (mineig <= 0) {
    Dmat <- Dmat + (abs(mineig) + 1e-6) * diag(m)
  }

  res <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq = 1, factorized = FALSE)
  w <- res$solution
  w[w < 0] <- 0
  w <- w / sum(w)

  return(list(w = w))
}

#' Constrained Least Squares Combiner (MC.Lse1+)
#'
#' Solves the unpenalized least squares regression combiner constrained to the probability simplex (\eqn{\lambda = 0}):
#' \deqn{\min_{\mathbf{w}} \|\mathbf{y} - \mathbf{F}\mathbf{w}\|_2^2 \quad \text{s.t.} \quad \sum_{m=1}^M w_m = 1, \; w_m \ge 0}
#'
#' @param y Numeric vector of responses of length \eqn{n}.
#' @param Fm Numeric matrix of base learner predictions of dimension \eqn{n \times M}.
#' @return A list containing \code{w} (estimated weight vector).
#' @export
lse_qp <- function(y, Fm) {
  ridge_qp(y, Fm, lambda = 0)
}

#' Accuracy-Based Weighted Average Combiner (MC.WAvg)
#'
#' Computes combiner weights inversely proportional to empirical prediction errors:
#' \deqn{w_m \propto 1 / \text{RMSE}_m}
#' normalized to sum to 1.
#'
#' @param y Numeric vector of responses of length \eqn{n}.
#' @param Fm Numeric matrix of base learner predictions of dimension \eqn{n \times M}.
#' @return A list containing \code{w} (normalized weight vector) and \code{rmse}.
#' @references
#' Wang, H., Fan, W., Yu, P. S., and Han, J. (2003). Mining concept-drifting data streams using ensemble classifiers.
#' In \emph{Proceedings of the Ninth ACM SIGKDD International Conference on Knowledge Discovery and Data Mining}, 226--235.
#' @export
weighted_avg <- function(y, Fm) {
  Fm <- as.matrix(Fm)
  y <- as.vector(y)
  M <- ncol(Fm)

  rmses <- numeric(M)
  for (m in seq_len(M)) {
    rmses[m] <- calc_rmse(Fm[, m], y)
  }

  inv_rmse <- 1 / (rmses + 1e-6)
  w <- inv_rmse / sum(inv_rmse)

  return(list(w = w, rmse = rmses))
}
