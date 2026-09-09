#' poisselect: Poisson Selection Models for Count Data
#'
#' Estimates a Heckman-type selection model for count data by maximum
#' likelihood. The whole estimation is implemented from scratch: the
#' Gauss-Hermite nodes and weights, the approximated log-likelihood on the log
#' scale together with its gradient, the optimisation with [stats::optim()] and
#' the standard errors from the numerical Hessian.
#'
#' @section Model:
#' The outcome equation models a count \eqn{y_i} that is observed for the
#' selected units only, the selection equation the fully observed indicator
#' \eqn{s_i}:
#' \deqn{\ln \mu_i = x_i' \beta + \varepsilon_i, \quad
#'       y_i \sim \mathrm{Poisson}(\mu_i),}
#' \deqn{s_i^* = z_i' \gamma + u_i, \quad s_i = 1\{s_i^* > 0\},}
#' with jointly normal error terms
#' \deqn{(\varepsilon_i, u_i) \sim N_2\left(0,
#'       \left[\begin{array}{cc} \sigma^2 & \rho\sigma \\
#'       \rho\sigma & 1 \end{array}\right]\right).}
#' The correlation \eqn{\rho} measures the strength of the selection bias; for
#' \eqn{\rho = 0} the selection is ignorable and a plain Poisson GLM on the
#' selected units would be consistent.
#'
#' @section Main entry points:
#' \describe{
#'   \item{[poisselect()]}{fits the model and returns an object of class
#'     `"poisselect"` with `print`, `summary`, `plot`, `predict`, `coef`,
#'     `vcov` and `logLik` methods.}
#'   \item{[simulate_poisselect()]}{simulates data from the model.}
#'   \item{[sim_selection]}{a ready-made simulated data set with known true
#'     parameters.}
#' }
#'
#' @keywords internal
#'
#' @importFrom checkmate assert_data_frame assert_formula assert_int
#'   assert_list assert_names assert_number assert_numeric assert_subset
#' @importFrom graphics abline barplot box legend par points
#' @importFrom stats .getXlevels binomial delete.response dnorm dpois glm.fit
#'   model.frame model.matrix model.response na.pass optim optimHess pnorm
#'   poisson printCoefmat quantile rnorm rpois setNames terms
#' @importFrom utils head modifyList
"_PACKAGE"
