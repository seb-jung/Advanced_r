# The main function only calls the named sub-tasks in order (top-down design).
# The argument 'K' keeps the upper-case symbol of the assignment; it is the
# only name that deviates from snake_case and is exempted from the linter.

#' Fit a Poisson Selection Model
#'
#' Estimates a Heckman-type selection model for count data by maximum
#' likelihood. A Poisson outcome equation with log-normal unobserved
#' heterogeneity is fitted jointly with a probit selection equation, so that
#' the correlation `rho` between the two error terms corrects the bias of a
#' plain Poisson GLM on the self-selected subsample.
#'
#' @details
#' The model consists of an outcome equation for the count `y`, which is
#' observed for the selected units only, and a selection equation for the fully
#' observed indicator `s`:
#' \deqn{\ln \mu_i = x_i'\beta + \varepsilon_i, \quad
#'       y_i \sim \mathrm{Poisson}(\mu_i),}
#' \deqn{s_i^* = z_i'\gamma + u_i, \quad s_i = 1\{s_i^* > 0\},}
#' \deqn{(\varepsilon_i, u_i) \sim N_2\left(0,
#'       \left[\begin{array}{cc} \sigma^2 & \rho\sigma \\
#'       \rho\sigma & 1\end{array}\right]\right).}
#'
#' The likelihood contains an integral over the outcome error that has no
#' closed form. It is approximated by a `K`-point Gauss-Hermite quadrature, and
#' the sum over the nodes is evaluated on the log scale with the log-sum-exp
#' trick, so that large counts and extreme selection probabilities cannot
#' overflow. The maximisation with [stats::optim()] runs on the unconstrained
#' parameters `log(sigma)` and `atanh(rho)` and uses the analytic gradient of
#' the log-likelihood. Standard errors come from the numerical Hessian at the
#' optimum and are transformed back to `sigma` and `rho` with the delta method.
#'
#' @section Accuracy of the quadrature:
#' A fixed `K`-point rule integrates accurately as long as the Poisson
#' probability is a smooth function of the outcome error, which is the case for
#' small and moderate counts. For counts in the hundreds combined with a large
#' `sigma`, `K = 20` can be too coarse and the approximated likelihood may even
#' have several local maxima. Refit with a larger `K` (say 40 or 80) and compare
#' the estimates and log-likelihoods in such cases.
#'
#' @section Identification:
#' The selection equation should contain at least one *exclusion restriction*,
#' that is a covariate which shifts the selection probability but does not
#' appear in the outcome equation. Without it the model is identified through
#' the functional form alone, which is weak in practice, and a warning is
#' issued.
#'
#' @param outcome A two-sided formula for the outcome equation, for example
#'   `y ~ x1 + x2`. The response must be a non-negative integer-valued count and
#'   should be `NA` for the non-selected units.
#' @param selection A two-sided formula for the selection equation, for example
#'   `s ~ x1 + z1`. The response must be a fully observed indicator with the
#'   values 0 and 1 (or a logical vector).
#' @param data A `data.frame` containing all variables of both formulas. The
#'   covariates of the selection equation must be complete for every unit, the
#'   covariates of the outcome equation only for the selected units.
#' @param K Number of Gauss-Hermite quadrature nodes, a single integer between
#'   2 and 200. The default of 20 is accurate for the usual range of counts, see
#'   the section on accuracy.
#' @param start Optional numeric vector of starting values on the original
#'   scale, ordered as `c(beta, gamma, sigma, rho)` with `sigma > 0` and
#'   `abs(rho) < 1`. Defaults to `NULL`, which derives the starting values from
#'   a probit fit of the selection equation and a Poisson GLM on the selected
#'   units.
#' @param control A named list of control parameters passed on to
#'   [stats::optim()], for example `list(maxit = 5000)`. Allowed entries are
#'   `maxit`, `reltol`, `abstol`, `trace`, `REPORT`, `fnscale`, `parscale` and
#'   `ndeps`.
#'
#' @return An object of class `"poisselect"`, a list with the components
#'   \describe{
#'     \item{`call`}{the matched call.}
#'     \item{`coefficients`}{a list with the outcome coefficients `outcome`
#'       (beta) and the selection coefficients `selection` (gamma).}
#'     \item{`sigma`, `rho`}{the estimated standard deviation of the outcome
#'       error and the estimated error correlation.}
#'     \item{`standard_errors`}{a list with the standard errors of `outcome`,
#'       `selection`, `sigma` and `rho`.}
#'     \item{`vcov`}{the covariance matrix of all parameters on the original
#'       scale.}
#'     \item{`loglik`, `npar`, `aic`}{the maximised log-likelihood, the number
#'       of estimated parameters and the AIC.}
#'     \item{`converged`, `convergence`, `counts`, `message`}{the convergence
#'       information reported by [stats::optim()].}
#'     \item{`n`, `n_selected`}{the number of units and the number of selected
#'       units.}
#'     \item{`K`, `quadrature`}{the number of quadrature nodes and the nodes
#'       and weights themselves.}
#'     \item{`model`}{the responses, design matrices, `terms` objects and
#'       factor levels of both equations, reused by the methods.}
#'     \item{`theta`}{the estimate on the unconstrained scale, reused by the
#'       profile plot.}
#'   }
#'
#' @seealso [summary.poisselect()], [plot.poisselect()],
#'   [predict.poisselect()], [simulate_poisselect()], [sim_selection]
#'
#' @examples
#' # A simulated data set with known parameters ships with the package.
#' fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)
#' fit
#'
#' # The true values are beta = c(0.5, 0.8, -0.4), gamma = c(0.3, 0.5, 0.7),
#' # sigma = 0.6 and rho = 0.5, see ?sim_selection.
#' summary(fit)
#'
#' # A plain Poisson GLM on the selected units is biased upwards here, because
#' # units with a large outcome error are also more likely to be selected.
#' selected <- subset(sim_selection, s == 1)
#' coef(glm(y ~ x1 + x2, family = poisson, data = selected))
#' coef(fit, which = "outcome")
#'
#' @export
poisselect <- function(outcome, selection, data,
                       K = 20, # nolint: object_name_linter.
                       start = NULL, control = list()) {
  check_arguments(outcome, selection, data, K, start, control)
  model <- build_model_data(outcome, selection, data)
  check_model_data(model)
  check_start_values(start, model)
  quadrature <- build_gauss_hermite(K)
  theta_start <- compute_start_values(model, start)
  optimum <- maximise_loglik(theta_start, model, quadrature, control)
  uncertainty <- compute_standard_errors(optimum$par, model, quadrature)
  new_poisselect(optimum, uncertainty, model, quadrature, match.call())
}
