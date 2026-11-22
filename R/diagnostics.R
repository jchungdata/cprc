#' Diagnostic and Metric Tools for CPRC
#'
#' Provides functions for concept drift quantification via angular divergence,
#' cluster extraction from combiner weights, and cumulative weight elbow diagnostics.
#'
#' @name diagnostics
NULL

#' Compute Concept Drift Angle Metric
#'
#' Computes the angular departure of an estimated combiner weight vector \eqn{\hat{\mathbf{w}}}
#' from the uniform weight vector \eqn{\mathbf{w}^* = (1/M, \dots, 1/M)^\top}:
#' \deqn{\eta(\mathbf{w}^*, \hat{\mathbf{w}}) = \arccos\left(\frac{|\langle \mathbf{w}^*, \hat{\mathbf{w}} \rangle|}{\|\mathbf{w}^*\|_2 \|\hat{\mathbf{w}}\|_2}\right)}
#'
#' @param w Numeric vector of estimated weights of length \eqn{M}. If a \code{"cprc_reg"} or
#'   \code{"cprc_class"} object is supplied, \code{w$w} is extracted.
#' @param M Integer. Number of base models (defaults to length of \code{w}).
#' @param unit Character. Either \code{"degrees"} (default) or \code{"radians"}.
#' @return Numeric scalar representing the concept drift angle.
#' @export
concept_drift_angle <- function(w, M = NULL, unit = c("degrees", "radians")) {
  unit <- match.arg(unit)
  if (is.list(w) && !is.null(w$w)) {
    w <- w$w
  }
  w <- as.vector(w)
  if (is.null(M)) M <- length(w)
  if (M < 2L) return(0)

  w_star <- rep(1 / M, M)
  norm_w <- sqrt(sum(w^2))
  norm_star <- sqrt(sum(w_star^2))

  if (norm_w <= 0 || norm_star <= 0) return(0)

  cos_val <- abs(sum(w * w_star)) / (norm_w * norm_star)
  cos_val <- pmin(pmax(cos_val, -1), 1)
  angle_rad <- acos(cos_val)

  if (unit == "degrees") {
    return(angle_rad * (180 / pi))
  } else {
    return(angle_rad)
  }
}

#' Extract Homogeneous Batch Clusters from CPRC Weights
#'
#' Groups base models that share approximately equal weights by thresholding pairwise differences.
#'
#' @param object A fitted model of class \code{"cprc_reg"} or \code{"cprc_class"}, or a numeric weight vector.
#' @param tol Numerical tolerance threshold for regarding two weights as equal. Default is \code{1e-3}.
#' @return A list containing \code{membership}, \code{csize}, \code{no}, and \code{unique_weights}.
#' @export
extract_clusters <- function(object, tol = 1e-3) {
  if (is.list(object) && !is.null(object$w)) {
    w <- object$w
  } else {
    w <- as.vector(object)
  }
  M <- length(w)
  if (M < 2L) {
    return(list(membership = 1L, csize = 1L, no = 1L, unique_weights = w))
  }

  adj <- matrix(0, M, M)
  for (i in seq_len(M - 1L)) {
    for (j in (i + 1L):M) {
      if (abs(w[i] - w[j]) <= tol) {
        adj[i, j] <- 1
        adj[j, i] <- 1
      }
    }
  }
  diag(adj) <- 1
  cl <- connected_components(adj)

  unique_w <- numeric(cl$no)
  for (g in seq_len(cl$no)) {
    members <- which(cl$membership == g)
    unique_w[g] <- mean(w[members])
  }

  return(list(
    membership = cl$membership,
    csize = cl$csize,
    no = cl$no,
    unique_weights = unique_w
  ))
}

#' Cumulative Weight Elbow Diagnostic Plot
#'
#' Generates a cumulative sum plot of sorted weights as introduced in Wang et al. (2018) to
#' detect sharp elbows that identify outlying subjects or distinct clusters.
#'
#' @param w Numeric vector of weights, or a fitted object of class \code{"cprc_unsupervised"},
#'   \code{"cprc_reg"}, or \code{"cprc_class"}.
#' @param sort_weights Logical. Whether to sort weights in ascending order. Default is \code{TRUE}.
#' @param transform Character. Either \code{"none"} (default) or \code{"exp"}.
#' @param main Title for the plot.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
#' @export
plot_elbow <- function(w, sort_weights = TRUE, transform = c("none", "exp"),
                       main = "Cumulative Weight Diagnostic Plot", ...) {
  transform <- match.arg(transform)
  if (is.list(w) && !is.null(w$w)) {
    w <- w$w
  }
  w <- as.vector(w)
  N <- length(w)

  w_vals <- if (sort_weights) sort(w) else w
  cum_w <- cumsum(w_vals) / sum(w_vals)

  if (transform == "exp") {
    y_vals <- exp(cum_w)
    ylab_text <- "exp(Cumulative Sum of Weights)"
  } else {
    y_vals <- cum_w
    ylab_text <- "Cumulative Sum of Weights"
  }

  graphics::plot(seq_len(N), y_vals, type = "b", pch = 19, col = "darkblue",
                 xlab = "Index of Sorted Weights", ylab = ylab_text,
                 main = main, ...)
  graphics::grid()
  invisible(list(indices = seq_len(N), values = y_vals))
}
