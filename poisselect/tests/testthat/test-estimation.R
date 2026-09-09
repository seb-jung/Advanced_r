# The most important tests of the package: if we simulate data from the
# model with known parameters, does poisselect() get those parameters back?
# The tolerances (0.12 for coefficients, 0.2 for rho) are chosen so that
# normal sampling noise at n = 2000 passes, but a sign error or a wrong
# scaling (e.g. a missing sqrt(2)) would fail.

test_that("the estimator recovers the parameters under selection bias", {
  fit <- fit_simulated(n = 2000L, rho = 0.5, seed = 11L)
  expect_true(fit$converged)
  expect_equal(unname(coef(fit, which = "outcome")), true_beta,
               tolerance = 0.12)
  expect_equal(unname(coef(fit, which = "selection")), true_gamma,
               tolerance = 0.12)
  expect_equal(fit$sigma, true_sigma, tolerance = 0.12)
  expect_equal(fit$rho, true_rho, tolerance = 0.2)
  expect_true(all(unlist(fit$standard_errors) > 0))
})

test_that("the estimator recovers the parameters without selection bias", {
  fit <- fit_simulated(n = 2000L, rho = 0, seed = 12L)
  expect_true(fit$converged)
  expect_equal(unname(coef(fit, which = "outcome")), true_beta,
               tolerance = 0.12)
  expect_equal(fit$sigma, true_sigma, tolerance = 0.12)
  # rho = 0 must not be rejected when the selection really is ignorable.
  expect_lt(abs(fit$rho / fit$standard_errors$rho), 2.5)
})

test_that("the estimator recovers a negative correlation", {
  fit <- fit_simulated(n = 2000L, rho = -0.6, seed = 13L)
  expect_true(fit$converged)
  expect_lt(fit$rho, 0)
  expect_equal(fit$rho, -0.6, tolerance = 0.2)
  expect_equal(unname(coef(fit, which = "outcome")), true_beta,
               tolerance = 0.12)
})

test_that("ignoring the selection biases a plain Poisson GLM", {
  # This is the whole point of the package. With rho > 0 the units with a
  # large outcome error are selected more often, so a plain GLM on the
  # selected rows overestimates the intercept. poisselect() should not.
  set.seed(21L)
  simulated <- simulate_poisselect(n = 2500L, beta = true_beta,
                                   gamma = true_gamma, sigma = true_sigma,
                                   rho = 0.7)
  fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = simulated)
  naive <- glm(y ~ x1 + x2, family = poisson, data = subset(simulated, s == 1))
  # Careful with the comparison: the GLM intercept also absorbs the
  # sigma^2 / 2 from the log-normal error, so we compare both on the scale
  # of log E[Y | x] = beta_0 + sigma^2 / 2 (see start_values.R).
  naive_intercept <- unname(coef(naive)[1L])
  corrected_intercept <- unname(coef(fit, which = "outcome")[1L]) +
    fit$sigma^2 / 2
  true_intercept <- true_beta[1L] + true_sigma^2 / 2
  expect_lt(abs(corrected_intercept - true_intercept),
            abs(naive_intercept - true_intercept))
  expect_gt(naive_intercept, true_intercept)
})

test_that("the Wald intervals cover the true parameters", {
  # Mini Monte Carlo: the 95 % Wald interval for rho should contain the true
  # rho in (almost) every replication. Six replications are not enough to
  # check the exact coverage, but enough to notice if the standard errors
  # are off by a factor, which is what a wrong delta method would do.
  covered <- vapply(seq_len(6L), function(replication) {
    fit <- fit_simulated(n = 1200L, rho = 0.5, seed = 100L + replication)
    interval <- fit$rho + c(-1, 1) * 1.96 * fit$standard_errors$rho
    true_rho >= interval[1L] && true_rho <= interval[2L]
  }, logical(1L))
  expect_gte(sum(covered), 5L)
})

test_that("the number of quadrature nodes barely changes the estimates", {
  coarse <- fit_simulated(n = 800L, rho = 0.5, seed = 31L, K = 20L)
  fine <- fit_simulated(n = 800L, rho = 0.5, seed = 31L, K = 60L)
  expect_equal(coef(coarse), coef(fine), tolerance = 0.01)
  # The change from K = 20 to K = 60 should be tiny compared to the standard
  # errors, otherwise K = 20 would not be a reasonable default.
  expect_true(all(abs(coef(coarse) - coef(fine)) <
                    0.1 * sqrt(diag(vcov(fine)))))
  expect_equal(coarse$loglik, fine$loglik, tolerance = 1e-4)
})

test_that("user-supplied starting values lead to the same optimum", {
  automatic <- fit_simulated(n = 800L, rho = 0.5, seed = 41L)
  manual <- fit_simulated(n = 800L, rho = 0.5, seed = 41L,
                          start = c(0.2, 0.5, -0.2, 0.1, 0.2, 0.4, 0.9, 0.1))
  expect_equal(coef(automatic), coef(manual), tolerance = 1e-3)
  expect_equal(automatic$loglik, manual$loglik, tolerance = 1e-6)
})

test_that("'control' is passed through to optim()", {
  # With maxit = 1 optim() cannot converge. If we get the non-convergence
  # warning, the control list really made it through to optim().
  expect_warning(
    fit <- fit_simulated(n = 400L, rho = 0.5, seed = 51L,
                         control = list(maxit = 1L)),
    "did not converge"
  )
  expect_false(fit$converged)
  expect_equal(fit$convergence, 1L)
})

test_that("the reported log-likelihood is a maximum of the likelihood", {
  fit <- reference_fit
  expect_equal(fit$loglik, compute_loglik(fit$theta, fit$model,
                                          fit$quadrature))
  # At a maximum the gradient should be (almost) zero, and moving any single
  # parameter a bit away should make the log-likelihood smaller.
  gradient <- compute_loglik_gradient(fit$theta, fit$model, fit$quadrature)
  expect_lt(max(abs(gradient)), 1e-3)
  perturbed <- vapply(seq_along(fit$theta), function(position) {
    theta <- fit$theta
    theta[position] <- theta[position] + 0.05
    compute_loglik(theta, fit$model, fit$quadrature)
  }, numeric(1L))
  expect_true(all(perturbed < fit$loglik))
})

test_that("boundary estimates and non-convergence are reported", {
  model <- build_model_data(y ~ x1 + x2, s ~ x1 + z1, sim_selection)
  at_rho_boundary <- list(convergence = 0L,
                          par = pack_parameters(true_beta, true_gamma,
                                                true_sigma, 0.999))
  expect_warning(warn_about_fit(at_rho_boundary, model), "'rho'.*boundary")
  at_sigma_boundary <- list(convergence = 0L,
                            par = pack_parameters(true_beta, true_gamma,
                                                  1e-5, true_rho))
  expect_warning(warn_about_fit(at_sigma_boundary, model), "'sigma'")
  failed <- list(convergence = 1L,
                 par = pack_parameters(true_beta, true_gamma, true_sigma,
                                       true_rho))
  expect_warning(warn_about_fit(failed, model), "did not converge")
})

test_that("a singular information matrix gives NA standard errors", {
  expect_warning(
    covariance <- invert_information(matrix(0, 2L, 2L)),
    "singular"
  )
  expect_true(all(is.na(covariance)))
  expect_warning(
    standard_errors <- extract_standard_errors(diag(c(1, -1))),
    "not positive definite"
  )
  expect_equal(standard_errors, c(1, NA))
})

test_that("the estimation works without an intercept", {
  set.seed(61L)
  simulated <- simulate_poisselect(n = 1500L, beta = c(0, 0.8, -0.4),
                                   gamma = c(0, 0.5, 0.7), rho = 0.4)
  fit <- poisselect(y ~ 0 + x1 + x2, s ~ 0 + x1 + z1, data = simulated)
  expect_equal(names(coef(fit, which = "outcome")), c("x1", "x2"))
  expect_equal(unname(coef(fit, which = "outcome")), c(0.8, -0.4),
               tolerance = 0.15)
})

test_that("the fit is reproducible and independent of the row order", {
  set.seed(71L)
  simulated <- simulate_poisselect(n = 500L, rho = 0.4)
  first <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = simulated, K = 10L)
  second <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = simulated, K = 10L)
  expect_identical(coef(first), coef(second))
  expect_identical(first$loglik, second$loglik)
  set.seed(72L)
  shuffled <- simulated[sample(nrow(simulated)), ]
  third <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = shuffled, K = 10L)
  expect_equal(coef(first), coef(third), tolerance = 1e-6)
  expect_equal(first$loglik, third$loglik, tolerance = 1e-8)
})
