#' Extract the Coefficients of a Poisson Selection Model
#'
#' @param object An object of class `"poisselect"`.
#' @param which Which coefficients to return: `"all"` (the default) returns the
#'   coefficients of both equations together with `sigma` and `rho`, using the
#'   prefixes `outcome_` and `selection_` to keep the names unambiguous;
#'   `"outcome"` and `"selection"` return the coefficients of one equation
#'   under their plain variable names.
#' @param ... Currently ignored, present for consistency with the generic.
#'
#' @return A named numeric vector.
#'
#' @examples
#' fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)
#' coef(fit)
#' coef(fit, which = "outcome")
#'
#' @export
coef.poisselect <- function(object, which = c("all", "outcome", "selection"),
                            ...) {
  which <- match.arg(which)
  switch(
    which,
    outcome = object$coefficients$outcome,
    selection = object$coefficients$selection,
    all = setNames(
      c(object$coefficients$outcome, object$coefficients$selection,
        object$sigma, object$rho),
      build_parameter_names(object$model)
    )
  )
}

#' Extract the Covariance Matrix of a Poisson Selection Model
#'
#' The matrix refers to the original parametrisation, that is to `sigma` and
#' `rho` rather than to `log(sigma)` and `atanh(rho)`, and covers all
#' parameters of both equations. Its row and column names match
#' `names(coef(object))`.
#'
#' @param object An object of class `"poisselect"`.
#' @param ... Currently ignored, present for consistency with the generic.
#'
#' @return A numeric matrix with one row and column per parameter.
#'
#' @examples
#' fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)
#' round(vcov(fit), 4)
#' # The standard errors are the square roots of the diagonal.
#' sqrt(diag(vcov(fit)))
#'
#' @export
vcov.poisselect <- function(object, ...) {
  object$vcov
}

#' Extract the Log-Likelihood of a Poisson Selection Model
#'
#' The returned object carries the `df` and `nobs` attributes, so
#' [stats::AIC()] and [stats::BIC()] work out of the box. Every unit
#' contributes to the likelihood, the non-selected ones through the selection
#' equation, so `nobs` is the full sample size.
#'
#' @param object An object of class `"poisselect"`.
#' @param ... Currently ignored, present for consistency with the generic.
#'
#' @return An object of class `"logLik"`.
#'
#' @examples
#' fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#'
#' @export
logLik.poisselect <- function(object, ...) {
  structure(object$loglik, df = object$npar, nobs = object$n,
            class = "logLik")
}
