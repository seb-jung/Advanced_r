#' Diagnostic Plots for a Poisson Selection Model
#'
#' **Plot 1 (fit of the count distribution)** compares the observed relative
#' frequencies of the counts among the selected units with the count
#' distribution that the fitted model implies for exactly those units. For a
#' single selected unit the implied distribution is
#' \deqn{\Pr(y_i = m \mid x_i, z_i, s_i = 1) \approx
#'   \frac{\sum_{k=1}^K \frac{w_k}{\sqrt\pi} p(m \mid \mu_{ik})
#'         \Phi(\eta_{ik})}
#'        {\sum_{k=1}^K \frac{w_k}{\sqrt\pi} \Phi(\eta_{ik})},}
#' and the plotted frequency of `m` is the average over all selected units.
#'
#' **Plot 2 (log-likelihood along rho)** evaluates the log-likelihood on a grid
#' of values for `rho` while all other parameters stay at their estimates. A
#' pronounced peak means that `rho` is well identified; a flat curve is the
#' typical symptom of a missing or weak exclusion restriction.
#'
#' @param x An object of class `"poisselect"`.
#' @param which Which plots to draw, a subset of `1:2`. Defaults to both.
#' @param max_count Largest count shown in the first plot. Defaults to `NULL`,
#'   which uses the largest observed count, capped at the 99th percentile so
#'   that a single extreme count cannot flatten the whole picture.
#' @param n_grid Number of grid points of the second plot, at least 5.
#' @param ... Currently ignored, present for consistency with the generic.
#'
#' @return `x`, invisibly. Called for the side effect of drawing.
#'
#' @examples
#' fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)
#' plot(fit)
#' plot(fit, which = 2)
#'
#' @export
plot.poisselect <- function(x, which = 1:2, max_count = NULL, n_grid = 41L,
                            ...) {
  assert_subset(which, choices = 1:2, empty.ok = FALSE, .var.name = "which")
  assert_int(max_count, lower = 1L, null.ok = TRUE, .var.name = "max_count")
  assert_int(n_grid, lower = 5L, upper = 1000L, .var.name = "n_grid")
  which <- unique(which)
  if (length(which) > 1L) {
    old_par <- par(mfrow = c(1L, 2L))
    on.exit(par(old_par), add = TRUE)
  }
  if (1L %in% which) {
    plot_count_distribution(x, max_count)
  }
  if (2L %in% which) {
    plot_rho_profile(x, n_grid)
  }
  invisible(x)
}

#' Plot 1: Observed versus Model-Implied Count Distribution
#'
#' @param object An object of class `"poisselect"`.
#' @param max_count Largest count to show, or `NULL`.
#'
#' @return `invisible(NULL)`.
#' @noRd
plot_count_distribution <- function(object, max_count) {
  counts <- compute_count_distribution(object, max_count)
  heights <- rbind(observed = counts$observed, implied = counts$implied)
  barplot(heights, beside = TRUE, names.arg = counts$count,
          col = c("grey70", "steelblue"), border = NA,
          xlab = "count", ylab = "relative frequency",
          main = "Observed vs. model-implied counts",
          ylim = c(0, max(heights) * 1.15))
  legend("topright", legend = c("observed", "model-implied"),
         fill = c("grey70", "steelblue"), border = NA, bty = "n")
  box()
  invisible(NULL)
}

#' Observed and Model-Implied Count Distribution of the Selected Units
#'
#' The denominator of the implied distribution does not depend on `m` and is
#' computed once; only the Poisson factor of the numerator changes with `m`.
#' Both are aggregated on the log scale like the likelihood itself.
#'
#' @param object An object of class `"poisselect"`.
#' @param max_count Largest count to consider, or `NULL` for the automatic
#'   choice.
#'
#' @return A list with the numeric vectors `count`, `observed` and `implied`.
#' @noRd
compute_count_distribution <- function(object, max_count) {
  observed_counts <- object$model$y_selected
  if (is.null(max_count)) {
    cap <- quantile(observed_counts, probs = 0.99, names = FALSE)
    max_count <- max(1L, min(max(observed_counts), ceiling(cap)))
  }
  support <- 0:max_count
  parameters <- list(beta = object$coefficients$outcome,
                     gamma = object$coefficients$selection,
                     sigma = object$sigma, rho = object$rho)
  pieces <- compute_node_pieces(parameters, object$model, object$quadrature)
  log_denominator <- log_sum_exp_rows(
    add_log_weight(pieces$log_selection, pieces$log_weight)
  )
  implied <- vapply(support, function(m) {
    log_numerator <- log_sum_exp_rows(compute_log_terms(m, pieces))
    mean(exp(log_numerator - log_denominator))
  }, numeric(1L))
  list(
    count = support,
    observed = tabulate(observed_counts + 1L, nbins = max_count + 1L) /
      length(observed_counts),
    implied = implied
  )
}

#' Plot 2: Log-Likelihood along rho
#'
#' @param object An object of class `"poisselect"`.
#' @param n_grid Number of grid points.
#'
#' @return `invisible(NULL)`.
#' @noRd
plot_rho_profile <- function(object, n_grid) {
  profile <- compute_rho_profile(object, n_grid)
  plot(profile$rho, profile$loglik, type = "l", lwd = 2, col = "steelblue",
       xlab = expression(rho), ylab = "log-likelihood",
       main = "Log-likelihood along rho")
  abline(v = object$rho, col = "firebrick", lwd = 2, lty = 2)
  abline(v = 0, col = "grey60", lty = 3)
  points(object$rho, object$loglik, pch = 19, col = "firebrick")
  legend("bottom", legend = c("log-likelihood", "estimate", "rho = 0"),
         col = c("steelblue", "firebrick", "grey60"), lty = 1:3,
         lwd = c(2, 2, 1), bty = "n")
  invisible(NULL)
}

#' Log-Likelihood on a Grid of Values for rho
#'
#' The grid stops short of the boundary, where `1 / sqrt(1 - rho^2)` diverges,
#' and always contains the estimate itself.
#'
#' @param object An object of class `"poisselect"`.
#' @param n_grid Number of grid points.
#'
#' @return A list with the numeric vectors `rho` and `loglik`.
#' @noRd
compute_rho_profile <- function(object, n_grid) {
  grid <- sort(unique(c(seq(-0.99, 0.99, length.out = n_grid), object$rho)))
  theta <- object$theta
  position <- length(theta)
  loglik <- vapply(grid, function(rho) {
    theta[position] <- atanh(rho)
    compute_loglik(theta, object$model, object$quadrature)
  }, numeric(1L))
  list(rho = grid, loglik = loglik)
}
