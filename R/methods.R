#' S3 Methods for cprc Package Objects
#'
#' Provides \code{print}, \code{summary}, \code{predict}, and \code{plot} methods
#' for \code{cprc_reg}, \code{cprc_class}, and \code{cprc_unsupervised} objects.
#'
#' @name methods
NULL

#' @export
print.cprc_reg <- function(x, ...) {
  cat("\nConstrained Penalized Regression Combiner (Regression)\n")
  cat("----------------------------------------------------\n")
  cat("Penalty:      ", x$penalty, "(", x$type, ")\n")
  cat("Tuning lambda:", format(round(x$lambda, 4)), "\n")
  cat("Drift angle:  ", format(round(x$angle, 2)), "degrees\n")
  cat("Models (M):   ", x$M, "\n")
  cat("Number of groups:", x$clusters$no, "\n")
  cat("\nEstimated Combiner Weights:\n")
  w_disp <- round(x$w, 4)
  names(w_disp) <- paste0("M", seq_along(w_disp))
  print(w_disp)
  cat("\n")
  invisible(x)
}

#' @export
summary.cprc_reg <- function(object, ...) {
  ans <- list(
    call = object$call,
    penalty = object$penalty,
    type = object$type,
    lambda = object$lambda,
    angle = object$angle,
    weights = object$w,
    clusters = object$clusters,
    rmse = calc_rmse(object$fitted_values, object$fitted_values + object$residuals),
    n = object$n,
    M = object$M
  )
  class(ans) <- "summary.cprc_reg"
  return(ans)
}

#' @export
print.summary.cprc_reg <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("CPRC Streaming Regression Summary\n")
  cat("---------------------------------\n")
  cat("Penalty:        ", x$penalty, "(", x$type, ")\n")
  cat("Lambda:         ", format(round(x$lambda, 4)), "\n")
  cat("Drift Angle:    ", format(round(x$angle, 2)), "degrees\n")
  cat("Training RMSE:  ", format(round(x$rmse, 4)), "\n")
  cat("Base Models (M):", x$M, "\n")
  cat("Unique Groups:  ", x$clusters$no, "\n\n")
  cat("Weight Distribution across Base Models:\n")
  w_df <- data.frame(
    Model = paste0("M", seq_along(x$weights)),
    Weight = round(x$weights, 4),
    Cluster = x$clusters$membership
  )
  print(w_df, row.names = FALSE)
  cat("\n")
  invisible(x)
}

#' @export
predict.cprc_reg <- function(object, newdata, ...) {
  if (missing(newdata)) {
    return(object$fitted_values)
  }
  newdata <- as.matrix(newdata)
  if (ncol(newdata) != length(object$w)) {
    stop("Number of columns in newdata must match length of combiner weight vector.")
  }
  drop(newdata %*% object$w)
}

#' @export
plot.cprc_reg <- function(x, type = c("weights", "residuals"), ...) {
  type <- match.arg(type)
  if (type == "weights") {
    M <- length(x$w)
    graphics::plot(seq_len(M), x$w, type = "h", lwd = 3, col = "steelblue",
                   xlab = "Base Model Index", ylab = "Weight",
                   main = paste0("CPRC Combiner Weights (", x$penalty, ", Angle = ", round(x$angle, 1), "\u00B0)"),
                   ylim = c(0, max(x$w) * 1.1), ...)
    graphics::points(seq_len(M), x$w, pch = 19, col = "darkblue")
    graphics::abline(h = 1 / M, col = "red", lty = 2)
  } else {
    graphics::plot(x$fitted_values, x$residuals, pch = 19, col = "darkgray",
                   xlab = "Fitted Values", ylab = "Residuals",
                   main = "Residuals vs Fitted Values", ...)
    graphics::abline(h = 0, col = "red", lty = 2)
  }
  invisible(x)
}

#' @export
print.cprc_class <- function(x, ...) {
  cat("\nConstrained Penalized Regression Combiner (Classification)\n")
  cat("---------------------------------------------------------\n")
  cat("Penalty:      ", x$penalty, "\n")
  cat("Tuning lambda:", format(round(x$lambda, 4)), "\n")
  cat("Drift angle:  ", format(round(x$angle, 2)), "degrees\n")
  cat("Models (M):   ", x$M, "\n")
  cat("Number of groups:", x$clusters$no, "\n")
  cat("\nEstimated Combiner Weights:\n")
  w_disp <- round(x$w, 4)
  names(w_disp) <- paste0("M", seq_along(w_disp))
  print(w_disp)
  cat("\n")
  invisible(x)
}

#' @export
summary.cprc_class <- function(object, ...) {
  ans <- list(
    call = object$call,
    penalty = object$penalty,
    lambda = object$lambda,
    angle = object$angle,
    weights = object$w,
    clusters = object$clusters,
    n = object$n,
    M = object$M
  )
  class(ans) <- "summary.cprc_class"
  return(ans)
}

#' @export
print.summary.cprc_class <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("CPRC Streaming Classification Summary\n")
  cat("--------------------------------------\n")
  cat("Penalty:        ", x$penalty, "\n")
  cat("Lambda:         ", format(round(x$lambda, 4)), "\n")
  cat("Drift Angle:    ", format(round(x$angle, 2)), "degrees\n")
  cat("Base Models (M):", x$M, "\n")
  cat("Unique Groups:  ", x$clusters$no, "\n\n")
  w_df <- data.frame(
    Model = paste0("M", seq_along(x$weights)),
    Weight = round(x$weights, 4),
    Cluster = x$clusters$membership
  )
  print(w_df, row.names = FALSE)
  cat("\n")
  invisible(x)
}

#' @export
predict.cprc_class <- function(object, newdata, type = c("class", "prob"), threshold = 0.5, ...) {
  type <- match.arg(type)
  if (missing(newdata)) {
    if (type == "prob") return(object$fitted_probabilities)
    return(object$fitted_classes)
  }
  newdata <- as.matrix(newdata)
  if (ncol(newdata) != length(object$w)) {
    stop("Number of columns in newdata must match length of combiner weight vector.")
  }
  h <- drop(newdata %*% object$w)
  prob <- 1 / (1 + exp(-pmax(pmin(h, 50), -50)))
  if (type == "prob") return(prob)
  return(ifelse(prob >= threshold, 1, 0))
}

#' @export
plot.cprc_class <- function(x, ...) {
  M <- length(x$w)
  graphics::plot(seq_len(M), x$w, type = "h", lwd = 3, col = "forestgreen",
                 xlab = "Base Classifier Index", ylab = "Weight",
                 main = paste0("CPRC Classification Combiner Weights (Angle = ", round(x$angle, 1), "\u00B0)"),
                 ylim = c(0, max(x$w) * 1.1), ...)
  graphics::points(seq_len(M), x$w, pch = 19, col = "darkgreen")
  graphics::abline(h = 1 / M, col = "red", lty = 2)
  invisible(x)
}

#' @export
print.cprc_unsupervised <- function(x, ...) {
  cat("\nRegularized Unsupervised Aggregation of Parameter Maps\n")
  cat("------------------------------------------------------\n")
  cat("Features/Voxels (M):", x$M, "\n")
  cat("Subjects/Sources (N):", x$N, "\n")
  cat("Tuning lambda:       ", format(round(x$lambda, 4)), "\n")
  cat("Iterations:          ", x$iterations, "\n")
  cat("Converged:           ", x$converged, "\n")
  cat("\nSubject Aggregation Weights:\n")
  w_disp <- round(x$w, 4)
  names(w_disp) <- paste0("Sub", seq_along(w_disp))
  print(w_disp)
  cat("\n")
  invisible(x)
}

#' @export
summary.cprc_unsupervised <- function(object, ...) {
  ans <- list(
    call = object$call,
    M = object$M,
    N = object$N,
    lambda = object$lambda,
    iterations = object$iterations,
    converged = object$converged,
    weights = object$w,
    weight_summary = summary(object$w)
  )
  class(ans) <- "summary.cprc_unsupervised"
  return(ans)
}

#' @export
print.summary.cprc_unsupervised <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Regularized Unsupervised Aggregation Summary\n")
  cat("--------------------------------------------\n")
  cat("Total Features/Voxels (M):", x$M, "\n")
  cat("Number of Subjects (N):   ", x$N, "\n")
  cat("Regularization Lambda:    ", format(round(x$lambda, 4)), "\n")
  cat("Converged:                ", x$converged, "in", x$iterations, "iterations\n\n")
  cat("Weight Quantiles across Subjects:\n")
  print(x$weight_summary)
  cat("\n")
  invisible(x)
}

#' @export
predict.cprc_unsupervised <- function(object, newdata, ...) {
  if (missing(newdata)) {
    return(object$combined_map)
  }
  newdata <- as.matrix(newdata)
  if (ncol(newdata) != length(object$w)) {
    stop("Number of columns in newdata must match number of subjects in object.")
  }
  drop(newdata %*% object$w)
}

#' @export
plot.cprc_unsupervised <- function(x, type = c("elbow", "weights", "objective"), ...) {
  type <- match.arg(type)
  if (type == "elbow") {
    plot_elbow(x$w, main = "Subject Cumulative Weight Elbow Plot", ...)
  } else if (type == "weights") {
    N <- length(x$w)
    graphics::plot(seq_len(N), x$w, type = "h", lwd = 3, col = "purple",
                   xlab = "Subject Index", ylab = "Weight",
                   main = "Subject Aggregation Weights",
                   ylim = c(0, max(x$w) * 1.1), ...)
    graphics::points(seq_len(N), x$w, pch = 19, col = "darkmagenta")
    graphics::abline(h = 1 / N, col = "red", lty = 2)
  } else if (type == "objective") {
    graphics::plot(seq_along(x$obj_history), x$obj_history, type = "b", pch = 19, col = "darkred",
                   xlab = "Iteration", ylab = "Objective Function Value",
                   main = "Optimization Objective History", ...)
  }
  invisible(x)
}
