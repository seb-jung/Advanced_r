# Shared fixtures. The reference fit on the bundled data is built once and
# reused by the method tests instead of being refitted in every file.

reference_fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)

# The true parameters behind sim_selection, see ?sim_selection.
true_beta <- c(0.5, 0.8, -0.4)
true_gamma <- c(0.3, 0.5, 0.7)
true_sigma <- 0.6
true_rho <- 0.5

# A small data set that really follows the model, so that every error a test
# observes is caused by the defect the test introduces.
make_valid_data <- function(n = 40L, seed = 99L) {
  set.seed(seed)
  simulated <- simulate_poisselect(n = n, rho = 0.3)
  simulated$y_complete <- NULL
  simulated
}

# Fit on freshly simulated data with the true parameters of sim_selection.
fit_simulated <- function(n, rho, sigma = true_sigma, seed = 1L, ...) {
  set.seed(seed)
  simulated <- simulate_poisselect(n = n, beta = true_beta, gamma = true_gamma,
                                   sigma = sigma, rho = rho)
  poisselect(y ~ x1 + x2, s ~ x1 + z1, data = simulated, ...)
}

# Central finite differences, used to check the analytic gradient.
numerical_gradient <- function(f, theta, step = 1e-5) {
  vapply(seq_along(theta), function(i) {
    up <- theta
    down <- theta
    up[i] <- up[i] + step
    down[i] <- down[i] - step
    (f(up) - f(down)) / (2 * step)
  }, numeric(1L))
}
