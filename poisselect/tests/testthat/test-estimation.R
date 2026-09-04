# The parameter recovery tests are the central validation: on data simulated
# from the model the estimator has to return the parameters that generated it.
# The tolerances survive ordinary sampling variation at these sample sizes but
# would catch a sign or scaling error.

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
  # With rho > 0 the units with a large outcome error are over-represented, so
  # the naive intercept is too large while poisselect() gets it right.
  set.seed(21L)
  simulated <- simulate_poisselect(n = 2500L, beta = true_beta,
                                   gamma = true_gamma, sigma = true_sigma,
                                   rho = 0.7)
  fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = simulated)
  naive <- glm(y ~ x1 + x2, family = poisson, data = subset(simulated, s == 1))
  # The naive intercept also absorbs exp(sigma^2 / 2), so the comparison is
  # made on the scale of the population mean.
  naive_intercept <- unname(coef(naive)[1L])
  corrected_intercept <- unname(coef(fit, which = "outcome")[1L]) +
    fit$sigma^2 / 2
  true_intercept <- true_beta[1L] + true_sigma^2 / 2
  expect_lt(abs(corrected_intercept - true_intercept),
            abs(naive_intercept - true_intercept))
  expect_gt(naive_intercept, true_intercept)
})

test_that("the Wald intervals cover the true parameters", {
  # Over repeated samples the 95 percent intervals should cover the truth
  # almost always; six replications catch standard errors that are off by a
  # factor, which is what a wrong delta method would cause.
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
  # The quadrature error must stay well below the statistical uncertainty.
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
  # One iteration cannot converge, which proves that the argument arrives.
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
  # The gradient vanishes at the optimum and every perturbation lowers the
  # log-likelihood.
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
