#' Construct an Object of Class "poisselect"
#'
#' @param optimum The list returned by [stats::optim()].
#' @param uncertainty The list returned by [compute_standard_errors()].
#' @param model A model list as built by [build_model_data()].
#' @param quadrature A list with the components `nodes` and `weights`.
#' @param call The matched call of [poisselect()].
#'
#' @return An object of class `"poisselect"`, see [poisselect()].
#' @noRd
new_poisselect <- function(optimum, uncertainty, model, quadrature, call) {
  parameters <- split_parameters(optimum$par, model$n_beta, model$n_gamma)
  names(parameters$beta) <- colnames(model$x)
  names(parameters$gamma) <- colnames(model$z)
  parameter_names <- build_parameter_names(model)
  covariance <- uncertainty$vcov
  dimnames(covariance) <- list(parameter_names, parameter_names)
  standard_errors <- uncertainty$standard_errors
  n_parameter <- length(parameter_names)
  # We minimised the negative log-likelihood, so optimum$value is -loglik.
  loglik <- -optimum$value
  object <- list(
    call = call,
    coefficients = list(outcome = parameters$beta,
                        selection = parameters$gamma),
    sigma = parameters$sigma,
    rho = parameters$rho,
    standard_errors = list(
      outcome = setNames(standard_errors[seq_len(model$n_beta)],
                         colnames(model$x)),
      selection = setNames(
        standard_errors[model$n_beta + seq_len(model$n_gamma)],
        colnames(model$z)
      ),
      sigma = standard_errors[n_parameter - 1L],
      rho = standard_errors[n_parameter]
    ),
    vcov = covariance,
    loglik = loglik,
    npar = n_parameter,
    aic = -2 * loglik + 2 * n_parameter,
    converged = optimum$convergence == 0L,
    convergence = optimum$convergence,
    counts = optimum$counts,
    message = optimum$message,
    n = model$n,
    n_selected = model$n_selected,
    K = length(quadrature$nodes),
    quadrature = quadrature,
    model = model,
    theta = optimum$par
  )
  class(object) <- "poisselect"
  object
}
