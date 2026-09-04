#' Starting Values for the Optimiser
#'
#' Both equations are first estimated separately and without any selection
#' correction, which is the model that would be right under `rho = 0`. The
#' starting value for `rho` is therefore 0, the hypothesis of ignorable
#' selection. A Poisson GLM estimates \eqn{\ln E[Y \mid x] = x'\beta +
#' \sigma^2/2}, so its intercept is shifted down by \eqn{\sigma_0^2/2} to
#' obtain a start for \eqn{\beta}.
#'
#' @param model A model list as built by [build_model_data()].
#' @param start User-supplied starting values on the original scale, ordered
#'   as `c(beta, gamma, sigma, rho)`, or `NULL`.
#'
#' @return Numeric vector of starting values on the unconstrained scale.
#' @noRd
compute_start_values <- function(model, start) {
  if (!is.null(start)) {
    n_coefficient <- model$n_beta + model$n_gamma
    return(pack_parameters(
      beta = start[seq_len(model$n_beta)],
      gamma = start[model$n_beta + seq_len(model$n_gamma)],
      sigma = start[n_coefficient + 1L],
      rho = start[n_coefficient + 2L]
    ))
  }
  gamma <- fit_probit_start(model)
  poisson_fit <- fit_glm_start(model$x_selected, model$y_selected, poisson())
  sigma <- estimate_start_sigma(model$y_selected, poisson_fit$fitted.values)
  beta <- poisson_fit$coefficients
  intercept <- match("(Intercept)", colnames(model$x), nomatch = 0L)
  if (intercept > 0L) {
    beta[intercept] <- beta[intercept] - sigma^2 / 2
  }
  pack_parameters(beta, gamma, sigma, 0)
}

#' Probit Starting Values with a Check for Perfect Separation
#'
#' If some combination of the selection covariates separates the selected from
#' the non-selected units perfectly, the probit coefficients have no finite
#' maximum likelihood estimate, and neither has `gamma` in the full model. The
#' probit fit reveals this through a deviance of numerically zero.
#'
#' @param model A model list as built by [build_model_data()].
#'
#' @return Numeric vector of length `n_gamma`.
#' @noRd
fit_probit_start <- function(model) {
  fit <- fit_glm_start(model$z, model$s, binomial(link = "probit"))
  if (fit$deviance < 1e-6 * model$n) {
    warning("The selection covariates separate the selected from the ",
            "non-selected units perfectly, so the coefficients of the ",
            "selection equation are not identified and their estimates and ",
            "standard errors are not meaningful.", call. = FALSE)
  }
  fit$coefficients
}

#' Fit a GLM for Starting Values Only
#'
#' The warnings of `glm.fit()` about fitted probabilities of 0 or 1 or slow
#' convergence are harmless for mere starting values and are silenced; any
#' non-finite coefficient is replaced by zero so that the optimiser never
#' receives NA.
#'
#' @param design A numeric design matrix.
#' @param response The corresponding response vector.
#' @param family A [stats::family()] object.
#'
#' @return The list returned by [stats::glm.fit()].
#' @noRd
fit_glm_start <- function(design, response, family) {
  fit <- suppressWarnings(glm.fit(x = design, y = response, family = family))
  fit$coefficients[!is.finite(fit$coefficients)] <- 0
  fit
}

#' Moment Estimator for the Starting Value of sigma
#'
#' Under the model \eqn{\mathrm{Var}(Y \mid x) = \mu + \mu^2 (e^{\sigma^2} -
#' 1)}, so the variation in excess of the Poisson benchmark identifies
#' \eqn{\sigma}. The value is bounded away from 0 and from implausibly large
#' values, because a degenerate starting value would stall the optimiser.
#'
#' @param y Numeric vector of observed counts.
#' @param fitted_mean Numeric vector of fitted Poisson means.
#'
#' @return A positive scalar.
#' @noRd
estimate_start_sigma <- function(y, fitted_mean) {
  excess <- mean((y - fitted_mean)^2 - fitted_mean) / mean(fitted_mean^2)
  sqrt(log1p(min(max(excess, 0.01), 10)))
}
