#' Simulate Data from a Poisson Selection Model
#'
#' Draws a data set that follows the model of [poisselect()] exactly, which
#' makes it possible to check whether the estimator recovers known parameters.
#'
#' @details
#' Three standard normal covariates are drawn: `x1` and `x2` enter the outcome
#' equation, `x1` and `z1` the selection equation, so `z1` is the exclusion
#' restriction that identifies the model. Both equations contain an intercept,
#' hence `beta` and `gamma` have three elements each, ordered as
#' `c(intercept, x1, x2)` and `c(intercept, x1, z1)`.
#'
#' The correlated errors are built from two independent standard normal draws
#' \eqn{u} and \eqn{e} as \eqn{\varepsilon = \sigma(\rho u +
#' \sqrt{1-\rho^2}\,e)}, which gives \eqn{\mathrm{Var}(\varepsilon) =
#' \sigma^2}, \eqn{\mathrm{Var}(u) = 1} and \eqn{\mathrm{Cov}(\varepsilon, u)
#' = \rho\sigma}. The function uses the random number generator of the current
#' session, so call [set.seed()] beforehand for reproducible data.
#'
#' @param n Number of units, a single integer of at least 2.
#' @param beta Numeric vector of length 3 with the coefficients of the outcome
#'   equation, ordered as `c(intercept, x1, x2)`.
#' @param gamma Numeric vector of length 3 with the coefficients of the
#'   selection equation, ordered as `c(intercept, x1, z1)`.
#' @param sigma Positive scalar, the standard deviation of the outcome error.
#' @param rho Scalar in `(-1, 1)`, the correlation of the two error terms.
#'
#' @return A `data.frame` with `n` rows and the columns
#'   \describe{
#'     \item{`y`}{the count outcome, `NA` for every non-selected unit.}
#'     \item{`s`}{the selection indicator, 0 or 1.}
#'     \item{`x1`, `x2`}{the covariates of the outcome equation.}
#'     \item{`z1`}{the exclusion restriction of the selection equation.}
#'     \item{`y_complete`}{the outcome of every unit, including the
#'       non-selected ones. Not used by the model, but handy for illustrating
#'       the size of the selection bias.}
#'   }
#'
#' @seealso [poisselect()], [sim_selection]
#'
#' @examples
#' set.seed(42)
#' simulated <- simulate_poisselect(n = 400, rho = 0.5)
#' str(simulated)
#' mean(simulated$s)
#'
#' fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = simulated)
#' coef(fit, which = "outcome")
#'
#' @export
simulate_poisselect <- function(n = 500L, beta = c(0.5, 0.8, -0.4),
                                gamma = c(0.3, 0.5, 0.7), sigma = 0.6,
                                rho = 0.5) {
  check_simulate_arguments(n, beta, gamma, sigma, rho)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  z1 <- rnorm(n)
  u <- rnorm(n)
  epsilon <- sigma * (rho * u + sqrt(1 - rho^2) * rnorm(n))
  s <- as.integer(gamma[1L] + gamma[2L] * x1 + gamma[3L] * z1 + u > 0)
  y_complete <- rpois(n, exp(beta[1L] + beta[2L] * x1 + beta[3L] * x2 +
                               epsilon))
  data.frame(
    y = ifelse(s == 1L, y_complete, NA_integer_),
    s = s,
    x1 = x1,
    x2 = x2,
    z1 = z1,
    y_complete = y_complete
  )
}

#' Check the Arguments of simulate_poisselect()
#'
#' @inheritParams simulate_poisselect
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_simulate_arguments <- function(n, beta, gamma, sigma, rho) {
  assert_int(n, lower = 2L, .var.name = "n")
  assert_numeric(beta, len = 3L, any.missing = FALSE, finite = TRUE,
                 .var.name = "beta")
  assert_numeric(gamma, len = 3L, any.missing = FALSE, finite = TRUE,
                 .var.name = "gamma")
  # checkmate has no open bounds, so sigma > 0 and |rho| < 1 are explicit.
  assert_number(sigma, finite = TRUE, .var.name = "sigma")
  if (sigma <= 0) {
    stop("'sigma' must be strictly positive, but is ", sigma, ".",
         call. = FALSE)
  }
  assert_number(rho, finite = TRUE, .var.name = "rho")
  if (abs(rho) >= 1) {
    stop("'rho' must lie strictly between -1 and 1, but is ", rho, ".",
         call. = FALSE)
  }
  invisible(TRUE)
}
