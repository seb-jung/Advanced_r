# Things that several test files need. testthat loads helper-*.R files
# before the tests, so everything defined here is available everywhere.
# The reference fit is computed once here instead of in every test file,
# which keeps the whole suite under 30 seconds.

reference_fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)

# The true parameters that were used to generate sim_selection (see
# ?sim_selection and data-raw/sim_selection.R).
true_beta <- c(0.5, 0.8, -0.4)
true_gamma <- c(0.3, 0.5, 0.7)
true_sigma <- 0.6
true_rho <- 0.5

# A small, valid data set for the input-check tests. We start from data that
# really follows the model and then break exactly one thing per test, so
# when an error appears we know it comes from that one defect.
make_valid_data <- function(n = 40L, seed = 99L) {
  set.seed(seed)
  simulated <- simulate_poisselect(n = n, rho = 0.3)
  simulated$y_complete <- NULL
  simulated
}

# Simulate a fresh data set with the true parameters from above and fit it.
# Extra arguments (K, start, control) go straight to poisselect().
fit_simulated <- function(n, rho, sigma = true_sigma, seed = 1L, ...) {
  set.seed(seed)
  simulated <- simulate_poisselect(n = n, beta = true_beta, gamma = true_gamma,
                                   sigma = sigma, rho = rho)
  poisselect(y ~ x1 + x2, s ~ x1 + z1, data = simulated, ...)
}

# Numerical gradient with central differences, (f(x + h) - f(x - h)) / 2h.
# Used to check that our analytic gradient is correct.
numerical_gradient <- function(f, theta, step = 1e-5) {
  vapply(seq_along(theta), function(i) {
    up <- theta
    down <- theta
    up[i] <- up[i] + step
    down[i] <- down[i] - step
    (f(up) - f(down)) / (2 * step)
  }, numeric(1L))
}
