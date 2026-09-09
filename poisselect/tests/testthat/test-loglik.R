# Our vectorised, log-scale likelihood is compared with a deliberately naive
# version that follows the formula from the assignment word by word: a loop
# over units, a loop over nodes, everything on the probability scale. The
# naive version is slow and would underflow for big counts, but for small
# test data it is the most convincing reference we have.

naive_loglik <- function(beta, gamma, sigma, rho, y, s, x, z, n_nodes) {
  quadrature <- build_gauss_hermite(n_nodes)
  contribution <- vapply(seq_along(s), function(i) {
    linear_selection <- sum(z[i, ] * gamma)
    if (s[i] == 0) {
      return(log(1 - pnorm(linear_selection)))
    }
    inner <- 0
    for (k in seq_len(n_nodes)) {
      rate <- exp(sum(x[i, ] * beta) + sqrt(2) * sigma * quadrature$nodes[k])
      eta <- (linear_selection + sqrt(2) * rho * quadrature$nodes[k]) /
        sqrt(1 - rho^2)
      inner <- inner + quadrature$weights[k] * dpois(y[i], rate) * pnorm(eta)
    }
    log(inner / sqrt(pi))
  }, numeric(1L))
  sum(contribution)
}

probit_loglik <- function(linear_selection, s) {
  sum(ifelse(s == 1,
             pnorm(linear_selection, log.p = TRUE),
             pnorm(linear_selection, lower.tail = FALSE, log.p = TRUE)))
}

make_fixture <- function(n = 60L, seed = 7L) {
  set.seed(seed)
  simulated <- simulate_poisselect(n = n, sigma = 0.4, rho = 0.3)
  build_model_data(y ~ x1 + x2, s ~ x1 + z1, simulated)
}

test_that("compute_loglik() matches the naive evaluation of the formula", {
  model <- make_fixture()
  beta <- c(0.4, 0.7, -0.3)
  gamma <- c(0.2, 0.4, 0.6)
  for (n_nodes in c(5L, 20L)) {
    expect_equal(
      compute_loglik(pack_parameters(beta, gamma, 0.5, 0.4), model,
                     build_gauss_hermite(n_nodes)),
      naive_loglik(beta, gamma, 0.5, 0.4, model$y, model$s, model$x, model$z,
                   n_nodes)
    )
  }
  # Sanity check of the test itself: if we feed the reference sigma / sqrt(2)
  # instead of sigma, the comparison must fail. Otherwise the test would
  # not notice a missing sqrt(2) in the real code either.
  expect_false(isTRUE(all.equal(
    compute_loglik(pack_parameters(beta, gamma, 0.5, 0.4), model,
                   build_gauss_hermite(20L)),
    naive_loglik(beta, gamma, 0.5 / sqrt(2), 0.4, model$y, model$s, model$x,
                 model$z, 20L)
  )))
})

test_that("the likelihood separates into two parts when rho = 0", {
  # Special case rho = 0: eta no longer depends on the node t_k, so Phi(eta)
  # can be pulled out of the sum and the likelihood splits into a plain
  # probit likelihood times a Poisson-log-normal mixture.
  model <- make_fixture()
  beta <- c(0.4, 0.7, -0.3)
  gamma <- c(0.2, 0.4, 0.6)
  quadrature <- build_gauss_hermite(20L)
  linear_outcome <- drop(model$x_selected %*% beta)
  mixture_loglik <- sum(vapply(seq_along(linear_outcome), function(i) {
    rate <- exp(linear_outcome[i] + sqrt(2) * 0.5 * quadrature$nodes)
    log(sum(quadrature$weights * dpois(model$y_selected[i], rate)) / sqrt(pi))
  }, numeric(1L)))
  expect_equal(
    compute_loglik(pack_parameters(beta, gamma, 0.5, 0), model, quadrature),
    probit_loglik(drop(model$z %*% gamma), model$s) + mixture_loglik
  )
})

test_that("the likelihood approaches the Poisson GLM as sigma goes to zero", {
  model <- make_fixture()
  beta <- c(0.4, 0.7, -0.3)
  gamma <- c(0.2, 0.4, 0.6)
  poisson_loglik <- sum(dpois(model$y_selected,
                              exp(drop(model$x_selected %*% beta)),
                              log = TRUE))
  expect_equal(
    compute_loglik(pack_parameters(beta, gamma, 1e-8, 0), model,
                   build_gauss_hermite(20L)),
    poisson_loglik + probit_loglik(drop(model$z %*% gamma), model$s),
    tolerance = 1e-6
  )
})

test_that("the likelihood stays finite for very large counts", {
  # y = 5000: dpois() is essentially 0 here, so without the log-sum-exp
  # trick we would get log(0) = -Inf.
  model <- make_fixture()
  model$y_selected[1L] <- 5000
  value <- compute_loglik(
    pack_parameters(c(0.4, 0.7, -0.3), c(0.2, 0.4, 0.6), 0.5, 0.4),
    model, build_gauss_hermite(20L)
  )
  expect_true(is.finite(value))

  # And a complete fit with counts in the hundreds has to run without any
  # warning, not just the single likelihood evaluation.
  set.seed(11L)
  large <- simulate_poisselect(n = 800L, beta = c(4.5, 0.4, -0.2),
                               sigma = 0.3, rho = 0.5)
  expect_gt(max(large$y, na.rm = TRUE), 150)
  fit <- expect_no_warning(poisselect(y ~ x1 + x2, s ~ x1 + z1, data = large))
  expect_true(is.finite(fit$loglik))
  expect_true(all(is.finite(unlist(fit$standard_errors))))
})

test_that("log_sum_exp_rows() survives magnitudes that would overflow", {
  a <- matrix(c(800, 799, -800, -801), nrow = 2L, byrow = TRUE)
  result <- log_sum_exp_rows(a)
  expect_equal(result, c(800, -800) + log1p(exp(-1)))
  small <- matrix(c(-1, -2, -3, 0.5, 1.5, -0.5), nrow = 2L)
  expect_equal(log_sum_exp_rows(small), log(rowSums(exp(small))))
})

test_that("compute_negative_loglik() flips the sign and caps overflow", {
  model <- make_fixture()
  quadrature <- build_gauss_hermite(20L)
  theta <- pack_parameters(c(0.4, 0.7, -0.3), c(0.2, 0.4, 0.6), 0.5, 0.4)
  expect_equal(compute_negative_loglik(theta, model, quadrature),
               -compute_loglik(theta, model, quadrature))
  # An absurd intercept (1e5) makes every Poisson probability 0 at every
  # node. We expect the penalty value 1e10, not Inf or NaN.
  broken <- pack_parameters(c(1e5, 0, 0), c(0.2, 0.4, 0.6), 0.5, 0.4)
  expect_equal(compute_negative_loglik(broken, model, quadrature), 1e10)
})

test_that("compute_node_pieces() returns the documented quantities", {
  model <- make_fixture()
  quadrature <- build_gauss_hermite(5L)
  parameters <- list(beta = c(0.4, 0.7, -0.3), gamma = c(0.2, 0.4, 0.6),
                     sigma = 0.5, rho = 0.4)
  pieces <- compute_node_pieces(parameters, model, quadrature)
  expect_equal(dim(pieces$log_rate), c(model$n_selected, 5L))
  expect_equal(dim(pieces$eta), c(model$n_selected, 5L))
  expect_equal(pieces$log_weight, log(quadrature$weights) - 0.5 * log(pi))
  # Check the formula log(mu_ik) = x_i'beta + sqrt(2) sigma t_k by hand for
  # the first unit.
  expect_equal(
    pieces$log_rate[1L, ],
    sum(model$x_selected[1L, ] * parameters$beta) +
      sqrt(2) * parameters$sigma * quadrature$nodes
  )
  # Same for eta_ik = (z_i'gamma + sqrt(2) rho t_k) / sqrt(1 - rho^2).
  expect_equal(
    pieces$eta[1L, ],
    (sum(model$z_selected[1L, ] * parameters$gamma) +
       sqrt(2) * parameters$rho * quadrature$nodes) /
      sqrt(1 - parameters$rho^2)
  )
  expect_equal(pieces$log_selection, pnorm(pieces$eta, log.p = TRUE))
})

test_that("the analytic gradient matches central differences", {
  model <- make_fixture(n = 120L, seed = 8L)
  quadrature <- build_gauss_hermite(12L)
  objective <- function(theta) compute_loglik(theta, model, quadrature)
  for (theta in list(
    pack_parameters(c(0.4, 0.7, -0.3), c(0.2, 0.4, 0.6), 0.5, 0.4),
    pack_parameters(c(0.1, 0.2, 0.1), c(-0.3, 0.9, -0.2), 1.3, -0.7)
  )) {
    expect_equal(compute_loglik_gradient(theta, model, quadrature),
                 numerical_gradient(objective, theta), tolerance = 1e-6)
  }
  expect_equal(compute_negative_gradient(theta, model, quadrature),
               -compute_loglik_gradient(theta, model, quadrature))
})

test_that("the parameter transformation is a bijection", {
  theta <- pack_parameters(c(1, 2), c(3, 4, 5), 0.7, -0.6)
  parameters <- split_parameters(theta, 2L, 3L)
  expect_equal(parameters$beta, c(1, 2))
  expect_equal(parameters$gamma, c(3, 4, 5))
  expect_equal(parameters$sigma, 0.7)
  expect_equal(parameters$rho, -0.6)
  expect_equal(
    pack_parameters(parameters$beta, parameters$gamma, parameters$sigma,
                    parameters$rho),
    theta
  )
})
