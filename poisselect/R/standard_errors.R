#' Covariance Matrix and Standard Errors at the Optimum
#'
#' The Hessian of the *negative* log-likelihood at the maximum is the observed
#' information, and its inverse is the asymptotic covariance matrix on the
#' unconstrained scale. Because the reported parameters are `sigma` and `rho`
#' rather than `log(sigma)` and `atanh(rho)`, that matrix is transformed with
#' the delta method: \eqn{V = J V_u J'} with the Jacobian \eqn{J} of the
#' back-transformation. \eqn{J} is diagonal with ones for the coefficients,
#' \eqn{d\sigma/d\ln\sigma = \sigma} and \eqn{d\rho/d\,\mathrm{atanh}\rho =
#' 1 - \rho^2}, so the product reduces to an elementwise scaling.
#'
#' @param theta Numeric vector of parameter estimates on the unconstrained
#'   scale.
#' @param model A model list as built by [build_model_data()].
#' @param quadrature A list with the components `nodes` and `weights`.
#'
#' @return A list with the covariance matrix `vcov` and the vector
#'   `standard_errors`, both on the original scale.
#' @noRd
compute_standard_errors <- function(theta, model, quadrature) {
  # optimHess() differentiates the analytic gradient numerically, which is
  # far more accurate than differencing the log-likelihood twice.
  information <- optimHess(
    theta,
    fn = compute_negative_loglik,
    gr = compute_negative_gradient,
    model = model,
    quadrature = quadrature
  )
  covariance <- invert_information(information)
  parameters <- split_parameters(theta, model$n_beta, model$n_gamma)
  jacobian <- c(rep(1, model$n_beta + model$n_gamma), parameters$sigma,
                1 - parameters$rho^2)
  covariance <- covariance * tcrossprod(jacobian)
  list(vcov = covariance, standard_errors = extract_standard_errors(covariance))
}

#' Invert the Observed Information Matrix
#'
#' @param information A square numeric matrix.
#'
#' @return The inverse, or a matrix of `NA` if the information is singular.
#' @noRd
invert_information <- function(information) {
  covariance <- tryCatch(solve(information), error = function(condition) NULL)
  if (is.null(covariance)) {
    warning("The Hessian at the optimum is singular, so no covariance matrix ",
            "could be computed and all standard errors are NA. This usually ",
            "points to a weakly identified model, for example a missing ",
            "exclusion restriction.", call. = FALSE)
    covariance <- matrix(NA_real_, nrow(information), ncol(information))
  }
  covariance
}

#' Standard Errors from a Covariance Matrix
#'
#' @param covariance A square numeric matrix.
#'
#' @return Numeric vector of standard errors, `NA` where the estimated variance
#'   is not positive.
#' @noRd
extract_standard_errors <- function(covariance) {
  variance <- diag(covariance)
  is_invalid <- !is.na(variance) & variance <= 0
  if (any(is_invalid)) {
    warning("The Hessian at the optimum is not positive definite, so ",
            sum(is_invalid), " standard error(s) are reported as NA. The ",
            "optimiser has most likely stopped at a saddle point or at a ",
            "boundary rather than at a maximum.", call. = FALSE)
    variance[is_invalid] <- NA_real_
  }
  sqrt(variance)
}
