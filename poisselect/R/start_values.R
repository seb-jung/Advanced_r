#' Starting Values for the Optimiser
#'
#' Good starting values matter a lot for BFGS. Our idea: fit both equations
#' separately with plain GLMs first (probit for the selection, Poisson on
#' the selected rows for the outcome). That is exactly the model under
#' rho = 0, so rho starts at 0 as well ("selection is ignorable" as the
#' starting hypothesis). glm()/glm.fit() is allowed for starting values
#' according to the assignment.
#'
#' One subtlety: with the log-normal error the Poisson GLM does not estimate
#' beta_0 but beta_0 + sigma^2 / 2 (because E[exp(eps)] = exp(sigma^2 / 2)).
#' So once we have a start for sigma we shift the intercept down by
#' sigma^2 / 2.
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
#' Perfect separation: if some combination of the z-variables splits the
#' selected and non-selected units perfectly (e.g. s = 1 exactly when
#' z1 > 0), the probit ML estimate does not exist, the coefficients run off
#' to infinity. The same then happens to gamma in the full model. The probit
#' fit is a cheap detector for this: its deviance goes to (numerically) 0.
#' Better to warn here than to return gamma = 408 with a standard error of
#' 10000.
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
#' glm.fit() likes to warn about "fitted probabilities of 0 or 1" or slow
#' convergence. For starting values that does not matter, so the warnings
#' are suppressed. If a coefficient comes back NA or Inf anyway we replace
#' it with 0, because optim() must never get a non-finite start.
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
#' Moment estimator for sigma from the overdispersion. Under the model
#' Var(Y | x) = mu + mu^2 * (exp(sigma^2) - 1), so whatever variance is left
#' after subtracting the Poisson part tells us something about sigma. We
#' clamp the result to a sensible range, because a start of sigma = 0 or
#' sigma = 5 would make optim() get stuck right away.
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
