#' Internal Utility Functions for cprc
#'
#' Helper functions for difference matrix construction, error metrics,
#' cross-validation fold partitioning, and graph connected components.
#'
#' @name utils
#' @keywords internal
NULL

#' Construct All-Pairs First-Order Difference Matrix
#'
#' Builds the \eqn{M' \times M} matrix \eqn{\mathbf{D}} with \eqn{M' = \binom{M}{2}}
#' where each row contains 1 at position \eqn{l} and -1 at position \eqn{m} (\eqn{l < m}).
#'
#' @param M Integer number of base models (\eqn{M \ge 2}).
#' @return Numeric matrix of dimension \eqn{\binom{M}{2} \times M}.
#' @keywords internal
diff_mat <- function(M) {
  if (M < 2) {
    stop("M must be at least 2 to construct a pairwise difference matrix.")
  }
  pairs <- utils::combn(M, 2)
  num_pair <- ncol(pairs)
  Dm <- matrix(0, nrow = num_pair, ncol = M)
  for (k in seq_len(num_pair)) {
    Dm[k, pairs[1, k]] <- 1
    Dm[k, pairs[2, k]] <- -1
  }
  return(Dm)
}

#' Root Mean Squared Error
#'
#' @param pred Numeric vector of predicted values.
#' @param obs Numeric vector of observed values.
#' @return Numeric scalar RMSE.
#' @keywords internal
calc_rmse <- function(pred, obs) {
  sqrt(mean((obs - pred)^2))
}

#' Binary Classification Misclassification Rate
#'
#' @param pred Numeric or integer vector of predicted class labels (0/1).
#' @param obs Numeric or integer vector of observed class labels (0/1).
#' @return Numeric scalar misclassification rate.
#' @keywords internal
calc_class_error <- function(pred, obs) {
  mean(pred != obs)
}

#' Generate K-Fold Cross-Validation Indices
#'
#' @param n Integer sample size.
#' @param K Integer number of folds.
#' @return Integer vector of length \eqn{n} assigning each observation to a fold in \code{1:K}.
#' @keywords internal
cv_folds <- function(n, K = 10) {
  K <- min(K, n)
  sample(rep(seq_len(K), length.out = n))
}

#' Connected Components of an Undirected Adjacency Matrix
#'
#' Computes connected components via breadth-first search without external graph dependencies.
#'
#' @param adj Symmetric binary adjacency matrix of dimension \eqn{M \times M}.
#' @return A list with \code{membership} (integer vector of cluster IDs),
#'   \code{csize} (integer vector of cluster sizes), and \code{no} (number of clusters).
#' @keywords internal
connected_components <- function(adj) {
  M <- nrow(adj)
  visited <- rep(FALSE, M)
  membership <- integer(M)
  csize <- integer(0)
  comp_id <- 0L

  for (i in seq_len(M)) {
    if (!visited[i]) {
      comp_id <- comp_id + 1L
      queue <- i
      visited[i] <- TRUE
      cnt <- 0L
      while (length(queue) > 0L) {
        curr <- queue[1L]
        queue <- queue[-1L]
        membership[curr] <- comp_id
        cnt <- cnt + 1L
        neighbors <- which(adj[curr, ] > 0 & !visited)
        if (length(neighbors) > 0L) {
          visited[neighbors] <- TRUE
          queue <- c(queue, neighbors)
        }
      }
      csize <- c(csize, cnt)
    }
  }

  return(list(membership = membership, csize = csize, no = comp_id))
}
