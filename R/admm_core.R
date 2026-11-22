#' Core Numerical Optimization Routines for CPRC via ADMM
#'
#' Implements proximal gradient descent, soft thresholding, exact projection onto
#' the probability simplex, and the Alternating Direction Method of Multipliers (ADMM).
#'
#' @name admm_core
NULL

#' Exact Projection onto the Probability Simplex
#'
#' Solves the Euclidean projection problem:
#' \deqn{\min_{\mathbf{w}} \frac{1}{2} \|\mathbf{w} - \mathbf{v}\|_2^2 \quad \text{s.t.} \quad \sum_{m=1}^M w_m = 1, \; w_m \ge 0}
#' using the fast \eqn{\mathcal{O}(M \log M)} sorting algorithm developed by Shalev-Shwartz and Singer (2006).
#'
#' @param v Numeric vector of length \eqn{M}.
#' @return Numeric vector \eqn{\mathbf{w}} of length \eqn{M} lying on the unit probability simplex.
#' @references
#' Shalev-Shwartz, S. and Singer, Y. (2006). Efficient learning of label ranking by soft projections onto polyhedra.
#' \emph{Journal of Machine Learning Research}, 7, 1567--1599.
#' @export
prox_simplex <- function(v) {
  nv <- length(v)
  if (nv == 1L) return(1)
  ord <- order(v, decreasing = TRUE)
  sv <- v[ord]
  csum <- cumsum(sv)
  tgam <- sv - (csum - 1) / seq_len(nv)

  pos_idx <- which(tgam > 0)
  if (length(pos_idx) == 0L) {
    pos <- 1L
  } else {
    pos <- max(pos_idx)
  }
  gamma <- (csum[pos] - 1) / pos

  w <- numeric(nv)
  w[ord[seq_len(pos)]] <- sv[seq_len(pos)] - gamma
  w[w < 0] <- 0
  sum_w <- sum(w)
  if (sum_w > 0) {
    w <- w / sum_w
  } else {
    w <- rep(1 / nv, nv)
  }
  return(w)
}

#' Loss Function Gradient
#'
#' Computes the gradient of the empirical loss function with respect to the combiner weights \eqn{\mathbf{w}}.
#'
#' @param tFm Transposed design matrix \eqn{\mathbf{X}^\top} of dimension \eqn{M \times n}.
#' @param y Response vector of length \eqn{n}.
#' @param h Linear predictor vector \eqn{\mathbf{h} = \mathbf{X}\mathbf{w}} of length \eqn{n}.
#' @param method Character. Either \code{"regression"} (squared error loss) or \code{"classification"} (negative log-likelihood).
#' @return Gradient vector of length \eqn{M}.
#' @keywords internal
loss_grad <- function(tFm, y, h, method = c("regression", "classification")) {
  method <- match.arg(method)
  if (method == "regression") {
    g <- drop(-tFm %*% (y - h))
  } else {
    prob <- 1 / (1 + exp(-pmax(pmin(h, 50), -50)))
    g <- drop(-tFm %*% (y - prob))
  }
  return(g)
}

#' Dual Optimization Update for L1 Norm
#'
#' Solves the elementwise dual subproblem under the \eqn{\ell_1} penalty:
#' \deqn{\min_{e_l} \lambda |e_l| + \frac{\rho}{2} (e_l - z_l)^2 \quad \text{s.t.} \quad |e_l| \le 1}
#'
#' @param zvec Numeric vector \eqn{\mathbf{z} = \mathbf{D}\mathbf{w}^{(t+1)} - \rho^{-1}\boldsymbol{\delta}^{(t)}}.
#' @param lambda Regularization parameter \eqn{\lambda \ge 0}.
#' @param rho Augmented Lagrangian penalty parameter \eqn{\rho > 0}.
#' @return Numeric vector \eqn{\mathbf{e}}.
#' @keywords internal
dual_optim_l1 <- function(zvec, lambda, rho) {
  evec <- sign(zvec) * pmax(abs(zvec) - lambda / rho, 0)
  flag <- abs(evec) >= 1
  evec[flag] <- sign(evec[flag])
  return(evec)
}

#' Dual Optimization Update for L2 Norm
#'
#' Solves the elementwise dual subproblem under the \eqn{\ell_2} penalty:
#' \deqn{\min_{e_l} \lambda e_l^2 + \frac{\rho}{2} (e_l - z_l)^2 \quad \text{s.t.} \quad |e_l| \le 1}
#'
#' @param zvec Numeric vector \eqn{\mathbf{z} = \mathbf{D}\mathbf{w}^{(t+1)} - \rho^{-1}\boldsymbol{\delta}^{(t)}}.
#' @param lambda Regularization parameter \eqn{\lambda \ge 0}.
#' @param rho Augmented Lagrangian penalty parameter \eqn{\rho > 0}.
#' @return Numeric vector \eqn{\mathbf{e}}.
#' @keywords internal
dual_optim_l2 <- function(zvec, lambda, rho) {
  evec <- (rho / (2 * lambda + rho)) * zvec
  flag <- abs(evec) >= 1
  evec[flag] <- sign(evec[flag])
  return(evec)
}

#' Proximal Gradient Step for Primal Update
#'
#' Minimizes the augmented Lagrangian with respect to \eqn{\mathbf{w} \in \mathcal{C}}:
#' \deqn{\min_{\mathbf{w} \in \mathcal{C}} l(\mathbf{w}) + \frac{\rho}{2} \|\mathbf{e}^{(t)} + \rho^{-1}\boldsymbol{\delta}^{(t)} - \mathbf{D}\mathbf{w}\|_2^2}
#' using proximal gradient linearization and projection onto the probability simplex.
#'
#' @param Fm Matrix \eqn{\mathbf{X}} of dimension \eqn{n \times M}.
#' @param tFm Matrix \eqn{\mathbf{X}^\top}.
#' @param y Vector of responses.
#' @param u Current weight vector \eqn{\mathbf{w}^{(t)}}.
#' @param rho Augmented Lagrangian parameter.
#' @param Dm2 Quadratic penalty matrix \eqn{\mathbf{D}^\top\mathbf{D}}.
#' @param vvec Auxiliary vector \eqn{\mathbf{D}^\top (\boldsymbol{\delta}/\rho + \mathbf{e})}.
#' @param L Lipschitz constant (maximum eigenvalue of Hessian bound).
#' @param method Loss type: \code{"regression"} or \code{"classification"}.
#' @param max_step Maximum inner iterations. Default is 100.
#' @param conv_eps Inner convergence tolerance. Default is 1e-10.
#' @return A list containing the updated weight vector \code{w}.
#' @keywords internal
pgm_step <- function(Fm, tFm, y, u, rho, Dm2, vvec, L, method, max_step = 100, conv_eps = 1e-10) {
  w <- u
  for (r in seq_len(max_step)) {
    hu <- drop(Fm %*% u)
    gu <- loss_grad(tFm, y, hu, method) + rho * drop(Dm2 %*% u - vvec)
    w <- prox_simplex(u - (1 / L) * gu)
    if (sum((w - u)^2) < conv_eps) {
      break
    }
    u <- w
  }
  return(list(w = w))
}

#' Core ADMM Algorithm for Constrained Penalized Regression Combiners
#'
#' Solves the constrained penalized regression combiner for a single regularization parameter \eqn{\lambda}:
#' \deqn{\min_{\mathbf{w}} l(\mathbf{w}) + \lambda \|\mathbf{D}\mathbf{w}\|_\gamma^\gamma \quad \text{s.t.} \quad \sum_{m=1}^M w_m = 1, \; w_m \ge 0}
#'
#' @param y Response vector of length \eqn{n}.
#' @param Fm Matrix of base-model predictions of dimension \eqn{n \times M}.
#' @param lambda Regularization tuning parameter \eqn{\lambda \ge 0}.
#' @param method Character. Either \code{"regression"} or \code{"classification"}.
#' @param norm Character. Penalty norm: \code{"l1"} (fused lasso) or \code{"l2"} (pairwise ridge).
#' @param Dm Penalty difference matrix. If \code{NULL}, generated via \code{diff_mat(M)}.
#' @param rho Augmented penalty parameter. Default is 100.
#' @param warm_start Optional list with starting values \code{list(w = NULL, delta = NULL, evec = NULL)}.
#' @param max_step Maximum ADMM iterations. Default is 5000.
#' @param conv_eps Convergence threshold. Default is 1e-8.
#' @param verbose Logical. Whether to print progress messages. Default is \code{FALSE}.
#' @return An object of class \code{"cprc_fit"} containing estimated weights, dual variables, and convergence status.
#' @export
admm_cprc <- function(y, Fm, lambda, method = c("regression", "classification"),
                      norm = c("l1", "l2"), Dm = NULL, rho = 100,
                      warm_start = list(w = NULL, delta = NULL, evec = NULL),
                      max_step = 5000, conv_eps = 1e-8, verbose = FALSE) {
  method <- match.arg(method)
  norm <- match.arg(norm)

  Fm <- as.matrix(Fm)
  y <- as.vector(y)
  n <- nrow(Fm)
  M <- ncol(Fm)

  if (M < 2L) {
    stop("Fm must contain at least 2 base models (columns).")
  }

  ABSTOL <- 1e-4
  RELTOL <- 1e-3

  dual_optim <- if (norm == "l1") dual_optim_l1 else dual_optim_l2
  tFm <- t(Fm)

  if (is.null(Dm)) {
    Dm <- diff_mat(M)
    num_pair <- nrow(Dm)
    tDm <- t(Dm)
    Dm2 <- M * diag(1, M, M) - matrix(1, M, M)
  } else {
    num_pair <- nrow(Dm)
    tDm <- t(Dm)
    Dm2 <- tDm %*% Dm
  }

  # Global Lipschitz constant L ensuring f(w) <= f_{L,u}(w)
  hess_scale <- if (method == "classification") 0.25 else 1.0
  L <- max(eigen(hess_scale * (tFm %*% Fm) + rho * Dm2, symmetric = TRUE, only.values = TRUE)$values)
  if (L <= 0) L <- 1

  # Initialize parameters
  w <- u <- if (is.null(warm_start$w)) rep(1 / M, M) else warm_start$w
  delta <- if (is.null(warm_start$delta)) rep(0, num_pair) else warm_start$delta
  evec <- if (is.null(warm_start$evec)) rep(0, num_pair) else warm_start$evec
  evec_prev <- evec

  conv <- FALSE
  for (iter in seq_len(max_step)) {
    # (a) Primal update via proximal gradient method
    vvec <- drop(tDm %*% (delta / rho + evec))
    sol <- pgm_step(Fm, tFm, y, u, rho, Dm2, vvec, L, method, max_step = 100, conv_eps = conv_eps)
    w <- sol$w

    # (b) Dual update
    diffw <- drop(Dm %*% w)
    zvec <- diffw - delta / rho
    evec <- dual_optim(zvec, lambda, rho)

    # (c) Lagrangian multiplier update
    delta <- -rho * (zvec - evec)

    # Primal and dual residual convergence checks (Boyd et al., 2011)
    r_norm <- sqrt(sum((diffw - evec)^2))
    s_norm <- sqrt(sum((rho * drop(tDm %*% (evec - evec_prev)))^2))

    eps_pri <- sqrt(num_pair) * ABSTOL + RELTOL * max(sqrt(sum(diffw^2)), sqrt(sum(evec^2)))
    eps_dual <- sqrt(M) * ABSTOL + RELTOL * sqrt(sum((drop(tDm %*% delta))^2))

    if (iter > 1L && r_norm < eps_pri && s_norm < eps_dual) {
      conv <- TRUE
      if (verbose) message(sprintf("ADMM converged at iteration %d (lambda = %.4f)", iter, lambda))
      break
    } else {
      evec_prev <- evec
      u <- w
    }
  }

  # Adjusted weights and batch cluster extraction
  # For l1 fused lasso, leverage exact soft-thresholded ADMM auxiliary variable evec;
  # for l2 pairwise ridge, threshold primal pairwise differences Dm %*% w.
  evec_final <- drop(Dm %*% w)
  cluster_diff <- if (norm == "l1") evec else evec_final
  tol_group <- if (norm == "l1") 1e-4 else (1 / M) * 1e-2
  adj_pairs <- utils::combn(M, 2)
  adj_mat <- matrix(0, M, M)

  for (k in seq_len(num_pair)) {
    if (abs(cluster_diff[k]) <= tol_group) {
      adj_mat[adj_pairs[1, k], adj_pairs[2, k]] <- 1
      adj_mat[adj_pairs[2, k], adj_pairs[1, k]] <- 1
    }
  }
  diag(adj_mat) <- 1
  cl <- connected_components(adj_mat)

  # Average weights within each homogeneous cluster
  w_adj <- w
  for (g in seq_len(cl$no)) {
    members <- which(cl$membership == g)
    w_adj[members] <- mean(w[members])
  }
  w_adj <- prox_simplex(w_adj)

  structure(
    list(
      w = w_adj,
      w_raw = w,
      clusters = cl,
      evec = evec_final,
      delta = delta,
      lambda = lambda,
      method = method,
      norm = norm,
      iter = iter,
      converged = conv
    ),
    class = "cprc_fit"
  )
}
