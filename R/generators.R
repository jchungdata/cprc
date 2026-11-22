#' Synthetic Data Stream Generators for Concept Drift Benchmarking
#'
#' Functions to generate synthetic streaming data batches with controlled gradual and abrupt
#' concept drift, noise, and recurrent concepts for testing ensemble combiners.
#'
#' @name generators
NULL

#' Generate Moving Hyperplane Streaming Data
#'
#' @param n1 Integer. Training sample size per batch. Default is 100.
#' @param n2 Integer. Test sample size per batch. Default is 500.
#' @param blocks Integer. Total number of batches. Default is 30.
#' @param flippoint Integer vector specifying batch indices where abrupt changes occur. Default is \code{c(15)}.
#' @param dimen Integer. Number of continuous features \eqn{p}. Default is 10.
#' @param delta Numeric vector of length 2: gradual and abrupt drift magnitudes. Default is \code{c(0.02, 0.2)}.
#' @param flipprob Numeric. Proportion of labels to randomly invert (noise). Default is 0.
#' @return A list with \code{train_batches}, \code{test_batches}, \code{change_points}, and \code{weights}.
#' @export
make_seq_plane <- function(n1 = 100, n2 = 500, blocks = 30, flippoint = c(15),
                           dimen = 10, delta = c(0.02, 0.2), flipprob = 0.0) {
  seq_train <- list()
  seq_test <- list()
  weights <- list()

  for (i in seq_len(blocks)) {
    x_tr <- matrix(stats::runif(n1 * dimen, 0, 1), nrow = n1, ncol = dimen)
    x_te <- matrix(stats::runif(n2 * dimen, 0, 1), nrow = n2, ncol = dimen)

    if (i == 1L) {
      w_i <- rep(1, dimen)
    } else if (i %in% flippoint) {
      w_i <- weights[[i - 1L]] + sample(c(-1, 1), dimen, replace = TRUE) * delta[2L]
    } else {
      w_i <- weights[[i - 1L]] + sample(c(-1, 1), dimen, replace = TRUE) * delta[1L]
    }
    weights[[i]] <- w_i
    threshold <- sum(w_i) / 2

    y_tr <- as.numeric(drop(x_tr %*% w_i) > threshold)
    y_te <- as.numeric(drop(x_te %*% w_i) > threshold)

    if (flipprob > 0) {
      flip_n <- as.integer(n1 * flipprob)
      if (flip_n > 0L) {
        flip_idx <- sample(seq_len(n1), flip_n)
        y_tr[flip_idx] <- 1 - y_tr[flip_idx]
      }
    }

    df_tr <- as.data.frame(x_tr)
    colnames(df_tr) <- paste0("X", seq_len(dimen))
    df_tr$y <- y_tr

    df_te <- as.data.frame(x_te)
    colnames(df_te) <- paste0("X", seq_len(dimen))
    df_te$y <- y_te

    seq_train[[i]] <- df_tr
    seq_test[[i]] <- df_te
  }

  return(list(
    train_batches = seq_train,
    test_batches = seq_test,
    change_points = flippoint,
    weights = weights
  ))
}

#' Generate Streaming Ensemble Algorithm (SEA) Benchmark Data
#'
#' @param batch_size Integer. Number of observations per batch. Default is 500.
#' @param num_batches Integer. Total number of batches. Default is 100.
#' @param theta Numeric vector of threshold values across concepts. Default is \code{c(8, 9, 7, 9.5)}.
#' @param noise_ratio Numeric. Proportion of class labels to randomly flip. Default is 0.10.
#' @return A list with \code{train_batches}, \code{test_batches}, \code{theta}, and \code{batches_per_concept}.
#' @export
make_seq_sea <- function(batch_size = 500, num_batches = 100,
                         theta = c(8, 9, 7, 9.5), noise_ratio = 0.10) {
  num_concepts <- length(theta)
  batches_per_concept <- ceiling(num_batches / num_concepts)

  train_list <- list()
  test_list <- list()

  for (b in seq_len(num_batches)) {
    concept_idx <- min(ceiling(b / batches_per_concept), num_concepts)
    th <- theta[concept_idx]

    x_tr <- matrix(stats::runif(batch_size * 3L, 0, 10), nrow = batch_size, ncol = 3L)
    y_tr <- ifelse(x_tr[, 1L] + x_tr[, 2L] <= th, 1, 0)
    if (noise_ratio > 0) {
      n_noise <- round(batch_size * noise_ratio)
      idx_noise <- sample(seq_len(batch_size), n_noise)
      y_tr[idx_noise] <- 1 - y_tr[idx_noise]
    }
    df_tr <- as.data.frame(x_tr)
    colnames(df_tr) <- c("X1", "X2", "X3")
    df_tr$y <- y_tr

    x_te <- matrix(stats::runif(1000L * 3L, 0, 10), nrow = 1000L, ncol = 3L)
    y_te <- ifelse(x_te[, 1L] + x_te[, 2L] <= th, 1, 0)
    df_te <- as.data.frame(x_te)
    colnames(df_te) <- c("X1", "X2", "X3")
    df_te$y <- y_te

    train_list[[b]] <- df_tr
    test_list[[b]] <- df_te
  }

  return(list(
    train_batches = train_list,
    test_batches = test_list,
    theta = theta,
    batches_per_concept = batches_per_concept
  ))
}

#' Generate Drifting Sequential Regression Data Stream
#'
#' @param n Integer. Sample size per batch. Default is 200.
#' @param M Integer. Total number of batches. Default is 20.
#' @param p Integer. Number of continuous predictors. Default is 10.
#' @param change_points Integer vector specifying batch indices where abrupt changes occur. Default is \code{c(10)}.
#' @param sigma Residual noise standard deviation. Default is 1.
#' @return A list with \code{seqdata}, \code{change_points}, and \code{alpha_history}.
#' @export
make_seq_reg <- function(n = 200, M = 20, p = 10, change_points = c(10), sigma = 1) {
  half_p <- floor(p / 2)
  alpha1 <- c(rep(5, half_p), rep(-2, p - half_p))
  alpha2 <- c(rep(-2, half_p), rep(5, p - half_p))

  seqdata <- list()
  alpha_history <- list()

  current_alpha <- alpha1
  for (m in seq_len(M)) {
    if (m %in% change_points) {
      current_alpha <- if (all(current_alpha == alpha1)) alpha2 else alpha1
    }
    alpha_history[[m]] <- current_alpha

    X <- matrix(stats::runif(n * p, 0, 1), nrow = n, ncol = p)
    y <- drop(X %*% current_alpha + stats::rnorm(n, 0, sigma))

    df <- as.data.frame(X)
    colnames(df) <- paste0("X", seq_len(p))
    df$y <- y
    seqdata[[m]] <- df
  }

  return(list(
    seqdata = seqdata,
    change_points = change_points,
    alpha_history = alpha_history
  ))
}
