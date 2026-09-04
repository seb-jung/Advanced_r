test_that("poisselect() returns a complete S3 object", {
  expect_s3_class(reference_fit, "poisselect")
  expect_true(all(
    c("call", "coefficients", "sigma", "rho", "standard_errors", "vcov",
      "loglik", "npar", "aic", "converged", "convergence", "n", "n_selected",
      "K", "quadrature", "model", "theta") %in% names(reference_fit)
  ))
  expect_equal(reference_fit$n, nrow(sim_selection))
  expect_equal(reference_fit$n_selected, sum(sim_selection$s == 1L))
  expect_equal(reference_fit$K, 20L)
  expect_equal(reference_fit$npar, 8L)
  expect_true(reference_fit$converged)
  expect_true(is.finite(reference_fit$loglik))
  expect_equal(reference_fit$call[[1L]], as.name("poisselect"))
  expect_gt(reference_fit$sigma, 0)
  expect_lt(abs(reference_fit$rho), 1)
})

test_that("coef() returns the requested coefficients", {
  all_parameters <- coef(reference_fit)
  expect_length(all_parameters, 8L)
  expect_equal(
    names(all_parameters),
    c("outcome_(Intercept)", "outcome_x1", "outcome_x2",
      "selection_(Intercept)", "selection_x1", "selection_z1", "sigma", "rho")
  )
  expect_equal(unname(all_parameters[7:8]),
               c(reference_fit$sigma, reference_fit$rho))
  outcome <- coef(reference_fit, which = "outcome")
  expect_equal(names(outcome), c("(Intercept)", "x1", "x2"))
  expect_equal(unname(outcome), unname(all_parameters[1:3]))
  selection <- coef(reference_fit, which = "selection")
  expect_equal(names(selection), c("(Intercept)", "x1", "z1"))
  expect_equal(unname(selection), unname(all_parameters[4:6]))
})

test_that("vcov() is a symmetric matrix matching the standard errors", {
  covariance <- vcov(reference_fit)
  expect_equal(dim(covariance), c(8L, 8L))
  expect_equal(dimnames(covariance),
               list(names(coef(reference_fit)), names(coef(reference_fit))))
  expect_equal(covariance, t(covariance))
  expect_true(all(diag(covariance) > 0))
  expect_equal(
    unname(sqrt(diag(covariance))),
    unname(unlist(reference_fit$standard_errors))
  )
})

test_that("the delta method scales sigma and rho correctly", {
  # The covariance matrix is reported on the original scale, so its entries
  # for sigma and rho are the unconstrained ones scaled by the derivatives
  # d sigma / d log(sigma) = sigma and d rho / d atanh(rho) = 1 - rho^2.
  unconstrained <- solve(optimHess(
    reference_fit$theta, compute_negative_loglik, compute_negative_gradient,
    model = reference_fit$model, quadrature = reference_fit$quadrature
  ))
  jacobian <- c(rep(1, 6L), reference_fit$sigma, 1 - reference_fit$rho^2)
  expect_equal(unname(vcov(reference_fit)),
               unconstrained * tcrossprod(jacobian))
})

test_that("logLik() carries df and nobs so that AIC and BIC work", {
  value <- logLik(reference_fit)
  expect_s3_class(value, "logLik")
  expect_equal(as.numeric(value), reference_fit$loglik)
  expect_equal(attr(value, "df"), 8L)
  expect_equal(attr(value, "nobs"), reference_fit$n)
  expect_equal(AIC(reference_fit), -2 * reference_fit$loglik + 2 * 8L)
  expect_equal(AIC(reference_fit), reference_fit$aic)
  expect_equal(BIC(reference_fit),
               -2 * reference_fit$loglik + log(reference_fit$n) * 8L)
})

test_that("print() reports the key numbers and returns invisibly", {
  output <- capture.output(print(reference_fit))
  for (pattern in c("Poisson selection model", "Observations: 800",
                    "selected: 474", "Log-likelihood", "AIC",
                    "Convergence: successful", "Outcome equation",
                    "Selection equation", "sigma", "rho")) {
    expect_true(any(grepl(pattern, output)), info = pattern)
  }
  capture.output(printed <- withVisible(print(reference_fit)))
  expect_false(printed$visible)
  expect_identical(printed$value, reference_fit)
})

test_that("summary() builds both coefficient tables with Wald tests", {
  summarised <- summary(reference_fit)
  expect_s3_class(summarised, "summary.poisselect")
  outcome_table <- summarised$coefficients$outcome
  expect_equal(dim(outcome_table), c(3L, 4L))
  expect_equal(colnames(outcome_table),
               c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(rownames(outcome_table), c("(Intercept)", "x1", "x2"))
  expect_equal(rownames(summarised$coefficients$selection),
               c("(Intercept)", "x1", "z1"))
  expect_equal(rownames(summarised$error_terms), c("sigma", "rho"))
  expect_equal(outcome_table[, "Estimate"], coef(reference_fit, "outcome"))
  expect_equal(unname(outcome_table[, "Std. Error"]),
               unname(reference_fit$standard_errors$outcome))
  expect_equal(outcome_table[, "z value"],
               outcome_table[, "Estimate"] / outcome_table[, "Std. Error"])
  expect_equal(outcome_table[, "Pr(>|z|)"],
               2 * pnorm(-abs(outcome_table[, "z value"])))
  expect_true(all(outcome_table[, "Pr(>|z|)"] >= 0 &
                    outcome_table[, "Pr(>|z|)"] <= 1))
  expect_equal(summarised$rho, reference_fit$rho)
  expect_equal(summarised$loglik, reference_fit$loglik)
  expect_equal(summarised$aic, AIC(reference_fit))
})

test_that("print.summary() shows rho, the log-likelihood and the AIC", {
  output <- capture.output(print(summary(reference_fit)))
  for (pattern in c("^rho = ", "AIC", "Log-likelihood", "Outcome equation",
                    "Selection equation", "Error terms", "Pr\\(>\\|z\\|\\)")) {
    expect_true(any(grepl(pattern, output)), info = pattern)
  }
  capture.output(printed <- withVisible(print(summary(reference_fit))))
  expect_false(printed$visible)
  expect_error(build_coefficient_table(c(a = 1), NA_real_), NA)
})

test_that("predict() returns the documented quantities", {
  link <- predict(reference_fit, type = "link")
  response <- predict(reference_fit, type = "response")
  pselect <- predict(reference_fit, type = "pselect")
  expect_length(link, reference_fit$n)
  expect_length(response, reference_fit$n)
  expect_length(pselect, reference_fit$n)
  expect_equal(unname(link),
               unname(drop(reference_fit$model$x %*%
                             coef(reference_fit, which = "outcome"))))
  # The unconditional population mean exp(x'beta + sigma^2 / 2).
  expect_equal(unname(response),
               unname(exp(link + reference_fit$sigma^2 / 2)))
  expect_true(all(response > exp(link)))
  expect_equal(unname(pselect),
               unname(pnorm(drop(reference_fit$model$z %*%
                                   coef(reference_fit, which = "selection")))))
  expect_true(all(pselect > 0 & pselect < 1))
  expect_equal(predict(reference_fit), response)
})

test_that("predict() honours newdata and keeps factor levels", {
  grid <- data.frame(x1 = c(-1, 0, 1), x2 = c(0, 0, 0), z1 = c(0, 0.5, 1))
  expect_length(predict(reference_fit, newdata = grid), 3L)
  expect_equal(
    unname(predict(reference_fit, newdata = grid, type = "link")),
    unname(drop(cbind(1, grid$x1, grid$x2) %*%
                  coef(reference_fit, which = "outcome")))
  )
  expect_equal(
    unname(predict(reference_fit, newdata = grid, type = "pselect")),
    unname(pnorm(drop(cbind(1, grid$x1, grid$z1) %*%
                        coef(reference_fit, which = "selection"))))
  )
  set.seed(81L)
  simulated <- simulate_poisselect(n = 600L, rho = 0.4)
  simulated$group <- factor(rep(c("a", "b"), length.out = nrow(simulated)))
  fit <- poisselect(y ~ x1 + group, s ~ x1 + z1, data = simulated, K = 10L)
  # newdata with a single level must still use the contrast coding of the fit.
  new_level <- data.frame(x1 = 0, group = factor("b", levels = c("a", "b")))
  expect_equal(unname(predict(fit, newdata = new_level, type = "link")),
               unname(sum(coef(fit, which = "outcome") * c(1, 0, 1))))
  expect_error(predict(fit, newdata = data.frame(x1 = 0, group = "zzz")),
               "new level")
})

test_that("plot() draws, returns invisibly and restores par()", {
  path <- tempfile(fileext = ".pdf")
  pdf(path)
  on.exit({
    dev.off()
    unlink(path)
  }, add = TRUE)
  before <- par("mfrow")
  expect_invisible(plot(reference_fit))
  expect_equal(par("mfrow"), before)
  expect_silent(plot(reference_fit, which = 1L))
  expect_silent(plot(reference_fit, which = 2L))
  expect_silent(plot(reference_fit, which = c(1L, 1L), max_count = 5L))
  expect_silent(plot(reference_fit, which = 2L, n_grid = 11L))
  expect_identical(plot(reference_fit, which = 1L), reference_fit)
})

test_that("the model-implied count distribution follows the formula", {
  counts <- compute_count_distribution(reference_fit, max_count = 40L)
  expect_equal(counts$count, 0:40)
  expect_true(all(counts$implied >= 0))
  # Truncating at 40 loses almost no mass.
  expect_gt(sum(counts$implied), 0.98)
  expect_lt(sum(counts$implied), 1 + 1e-8)
  observed_counts <- sim_selection$y[sim_selection$s == 1L]
  expect_equal(counts$observed[1L], mean(observed_counts == 0L))
  expect_lt(max(abs(counts$observed - counts$implied)), 0.05)

  # Naive evaluation of the assignment's formula for m = 2, averaged over the
  # selected units.
  quadrature <- reference_fit$quadrature
  model <- reference_fit$model
  beta <- coef(reference_fit, "outcome")
  gamma <- coef(reference_fit, "selection")
  naive <- mean(vapply(seq_len(model$n_selected), function(i) {
    rate <- exp(sum(model$x_selected[i, ] * beta) +
                  sqrt(2) * reference_fit$sigma * quadrature$nodes)
    eta <- (sum(model$z_selected[i, ] * gamma) +
              sqrt(2) * reference_fit$rho * quadrature$nodes) /
      sqrt(1 - reference_fit$rho^2)
    weight <- quadrature$weights / sqrt(pi) * pnorm(eta)
    sum(weight * dpois(2, rate)) / sum(weight)
  }, numeric(1L)))
  expect_equal(counts$implied[3L], naive)

  # Without max_count a single extreme count does not stretch the support.
  stretched <- reference_fit
  stretched$model$y_selected[1L] <- 500L
  expect_lt(length(compute_count_distribution(stretched, NULL)$count), 100L)
})

test_that("the rho profile peaks at the estimate", {
  profile <- compute_rho_profile(reference_fit, 41L)
  expect_true(all(is.finite(profile$loglik)))
  expect_false(is.unsorted(profile$rho))
  expect_true(reference_fit$rho %in% profile$rho)
  expect_equal(profile$rho[which.max(profile$loglik)], reference_fit$rho)
  expect_equal(max(profile$loglik), reference_fit$loglik)
})
