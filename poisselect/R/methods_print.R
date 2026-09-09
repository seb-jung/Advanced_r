#' Print a Poisson Selection Model
#'
#' Shows the call, the sample sizes, the maximised log-likelihood, the
#' convergence status and the estimated coefficients of both equations together
#' with `sigma` and `rho`.
#'
#' @param x An object of class `"poisselect"`.
#' @param digits Number of significant digits for the printed numbers.
#' @param ... Currently ignored, present for consistency with the generic.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)
#' print(fit)
#'
#' @export
print.poisselect <- function(x, digits = max(3L, getOption("digits") - 3L),
                             ...) {
  cat("\nPoisson selection model\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nObservations: ", x$n, " (selected: ", x$n_selected,
      ", not selected: ", x$n - x$n_selected, ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = digits), " on ", x$npar,
      " parameters   AIC: ", format(x$aic, digits = digits), "\n", sep = "")
  cat("Quadrature nodes: ", x$K, "\n", sep = "")
  cat("Convergence: ", describe_convergence(x), "\n", sep = "")
  print_coefficient_block("Outcome equation (beta)", x$coefficients$outcome,
                          digits)
  print_coefficient_block("Selection equation (gamma)",
                          x$coefficients$selection, digits)
  print_coefficient_block("Error terms", c(sigma = x$sigma, rho = x$rho),
                          digits)
  cat("\n")
  invisible(x)
}

#' Print One Block of Named Coefficients
#'
#' @param heading Character string printed above the block.
#' @param coefficient Named numeric vector.
#' @param digits Number of significant digits.
#'
#' @return `invisible(NULL)`.
#' @noRd
print_coefficient_block <- function(heading, coefficient, digits) {
  cat("\n", heading, ":\n", sep = "")
  print.default(format(coefficient, digits = digits), print.gap = 2L,
                quote = FALSE)
  invisible(NULL)
}

#' Describe the Convergence Status in Words
#'
#' @param object An object of class `"poisselect"` or `"summary.poisselect"`.
#'
#' @return A single character string.
#' @noRd
describe_convergence <- function(object) {
  if (object$converged) {
    return("successful (optim code 0)")
  }
  paste0("FAILED (optim code ", object$convergence, ")")
}
