#' Maximise the Approximated Log-Likelihood
#'
#' Runs BFGS on the unconstrained scale with the analytic gradient. A tight
#' relative tolerance is affordable because every evaluation is vectorised, and
#' it matters because the Hessian is taken at the returned point.
#'
#' @param theta_start Numeric vector of starting values on the unconstrained
#'   scale.
#' @param model A model list as built by [build_model_data()].
#' @param quadrature A list with the components `nodes` and `weights`.
#' @param control A named list of control parameters for [stats::optim()].
#'
#' @return The list returned by [stats::optim()].
#' @noRd
maximise_loglik <- function(theta_start, model, quadrature, control) {
  control <- modifyList(list(maxit = 1000L, reltol = 1e-10), control)
  optimum <- optim(
    par = theta_start,
    fn = compute_negative_loglik,
    gr = compute_negative_gradient,
    method = "BFGS",
    control = control,
    model = model,
    quadrature = quadrature
  )
  warn_about_fit(optimum, model)
  optimum
}

#' Warn About Non-Convergence and Boundary Estimates
#'
#' Both are warnings rather than errors: the returned object is still useful
#' for diagnosis, the user just must not trust it blindly.
#'
#' @param optimum The list returned by [stats::optim()].
#' @param model A model list as built by [build_model_data()].
#'
#' @return `invisible(TRUE)`.
#' @noRd
warn_about_fit <- function(optimum, model) {
  if (optimum$convergence != 0L) {
    warning("optim() did not converge (code ", optimum$convergence, "). The ",
            "estimates should not be trusted. Try a larger 'maxit' in ",
            "'control' or other starting values via 'start'.", call. = FALSE)
  }
  parameters <- split_parameters(optimum$par, model$n_beta, model$n_gamma)
  if (abs(parameters$rho) > 0.99) {
    warning("The estimate of 'rho' lies at the boundary of the parameter ",
            "space (|rho| > 0.99); the standard errors are not reliable.",
            call. = FALSE)
  }
  if (parameters$sigma < 1e-3) {
    warning("The estimate of 'sigma' is practically zero, i.e. at the ",
            "boundary of the parameter space; the standard errors are not ",
            "reliable.", call. = FALSE)
  }
  invisible(TRUE)
}
