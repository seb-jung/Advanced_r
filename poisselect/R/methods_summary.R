#' Summarise a Poisson Selection Model
#'
#' Builds the full coefficient tables of both equations. Each row reports the
#' estimate, its standard error, the Wald statistic `z = Estimate / SE` and the
#' two-sided p-value of the hypothesis that the coefficient is zero, using the
#' asymptotic standard normal distribution of `z` under that hypothesis. The
#' error parameters `sigma` and `rho` get the same treatment; the row for `rho`
#' is a direct test for selection bias, because `rho = 0` means that the
#' selection is ignorable.
#'
#' @param object An object of class `"poisselect"`.
#' @param ... Currently ignored, present for consistency with the generic.
#'
#' @return An object of class `"summary.poisselect"`, a list with the
#'   coefficient tables `coefficients$outcome`, `coefficients$selection` and
#'   `error_terms`, plus `sigma`, `rho`, the log-likelihood, the AIC, the
#'   sample sizes and the convergence information.
#'
#' @examples
#' fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)
#' summary(fit)
#'
#' # The coefficient tables are ordinary matrices.
#' summary(fit)$coefficients$outcome
#'
#' @export
summary.poisselect <- function(object, ...) {
  result <- list(
    call = object$call,
    coefficients = list(
      outcome = build_coefficient_table(object$coefficients$outcome,
                                        object$standard_errors$outcome),
      selection = build_coefficient_table(object$coefficients$selection,
                                          object$standard_errors$selection)
    ),
    error_terms = build_coefficient_table(
      c(sigma = object$sigma, rho = object$rho),
      c(object$standard_errors$sigma, object$standard_errors$rho)
    ),
    sigma = object$sigma,
    rho = object$rho,
    loglik = object$loglik,
    aic = object$aic,
    npar = object$npar,
    n = object$n,
    n_selected = object$n_selected,
    K = object$K,
    converged = object$converged,
    convergence = object$convergence
  )
  class(result) <- "summary.poisselect"
  result
}

#' Build a Coefficient Table with Wald Tests
#'
#' @param estimate Named numeric vector of point estimates.
#' @param standard_error Numeric vector of standard errors.
#'
#' @return Numeric matrix with the columns `Estimate`, `Std. Error`, `z value`
#'   and `Pr(>|z|)`.
#' @noRd
build_coefficient_table <- function(estimate, standard_error) {
  z_value <- estimate / standard_error
  # Two-sided p-value 2 * (1 - Phi(|z|)), evaluated as the upper tail.
  p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
  table <- cbind(estimate, standard_error, z_value, p_value)
  dimnames(table) <- list(names(estimate),
                          c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  table
}

#' Print the Summary of a Poisson Selection Model
#'
#' @param x An object of class `"summary.poisselect"`.
#' @param digits Number of significant digits for the printed numbers.
#' @param signif_stars Logical flag, whether to mark the p-values with
#'   significance stars. Defaults to the `show.signif.stars` option.
#' @param ... Passed on to [stats::printCoefmat()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)
#' print(summary(fit), digits = 3)
#'
#' @export
print.summary.poisselect <- function(x,
                                     digits = max(3L, getOption("digits") - 3L),
                                     signif_stars = getOption(
                                       "show.signif.stars"
                                     ),
                                     ...) {
  cat("\nPoisson selection model\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nObservations: ", x$n, " (selected: ", x$n_selected,
      ", not selected: ", x$n - x$n_selected, ")\n", sep = "")
  print_coefficient_table("Outcome equation (beta)", x$coefficients$outcome,
                          digits, signif_stars, ...)
  print_coefficient_table("Selection equation (gamma)",
                          x$coefficients$selection, digits, signif_stars, ...)
  print_coefficient_table("Error terms", x$error_terms, digits, signif_stars,
                          ...)
  cat("\nrho = ", format(x$rho, digits = digits), " (",
      describe_selection_bias(x$rho), ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = digits), " on ", x$npar,
      " parameters   AIC: ", format(x$aic, digits = digits), "\n", sep = "")
  cat("Quadrature nodes: ", x$K, "\n", sep = "")
  cat("Convergence: ", describe_convergence(x), "\n\n", sep = "")
  invisible(x)
}

#' Print One Coefficient Table
#'
#' @param heading Character string printed above the table.
#' @param table A coefficient matrix as built by [build_coefficient_table()].
#' @param digits Number of significant digits.
#' @param signif_stars Logical flag, whether to print significance stars.
#' @param ... Passed on to [stats::printCoefmat()].
#'
#' @return `invisible(NULL)`.
#' @noRd
print_coefficient_table <- function(heading, table, digits, signif_stars, ...) {
  cat("\n", heading, ":\n", sep = "")
  printCoefmat(table, digits = digits, signif.stars = signif_stars,
               has.Pvalue = TRUE, P.values = TRUE, na.print = "NA", ...)
  invisible(NULL)
}

#' Describe the Direction of the Selection Bias
#'
#' @param rho The estimated error correlation.
#'
#' @return A single character string.
#' @noRd
describe_selection_bias <- function(rho) {
  if (abs(rho) < 1e-4) {
    return("selection is effectively ignorable")
  }
  if (rho > 0) {
    return("units with a large outcome are selected more often")
  }
  "units with a large outcome are selected less often"
}
