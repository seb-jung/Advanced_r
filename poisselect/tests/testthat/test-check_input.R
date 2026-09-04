# Every check gets its own expectation, and every expectation pins down a part
# of the message, so a check cannot silently be replaced by a different one.

test_that("the formulas must be two-sided formulas without '.'", {
  valid <- make_valid_data()
  expect_error(poisselect("y ~ x1", s ~ x1 + z1, valid), "outcome.*formula")
  expect_error(poisselect(y ~ x1, "s ~ z1", valid), "selection.*formula")
  expect_error(poisselect(NULL, s ~ x1 + z1, valid), "outcome")
  expect_error(poisselect(~ x1 + x2, s ~ x1 + z1, valid),
               "'outcome' must be a two-sided formula")
  expect_error(poisselect(y ~ x1 + x2, ~ x1 + z1, valid),
               "'selection' must be a two-sided formula")
  expect_error(poisselect(y ~ ., s ~ x1 + z1, valid), "'\\.' shortcut")
  expect_error(poisselect(y ~ x1, s ~ ., valid), "'\\.' shortcut")
})

test_that("'data' must be a data.frame that contains all variables", {
  valid <- make_valid_data()
  expect_error(poisselect(y ~ x1, s ~ x1 + z1), "'data' is missing")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, data = as.matrix(valid)),
               "data.frame")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, data = as.list(valid)),
               "data.frame")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, data = NULL), "data.frame")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, data = valid[0L, ]), "rows")
  expect_error(poisselect(y ~ x1 + nonsense, s ~ x1 + z1, valid),
               "not columns of 'data'.*'nonsense'")
  expect_error(poisselect(y ~ x1, s ~ x1 + missing_z, valid), "'missing_z'")
})

test_that("'K' must be a single integer between 2 and 200", {
  valid <- make_valid_data()
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, K = 1L), "K.*>= 2")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, K = 0L), "K")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, K = -5L), "K")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, K = 2.5), "K")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, K = c(10L, 20L)), "K")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, K = "20"), "K")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, K = NA), "K")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, K = 5000L), "K.*<= 200")
})

test_that("'control' must be a named list of optim() controls", {
  valid <- make_valid_data()
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, control = "maxit"),
               "control")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, control = list(1000)),
               "names")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid,
                          control = list(bogus = 1)),
               "subset")
})

test_that("'start' must have the right length and stay in range", {
  valid <- make_valid_data()
  # 2 outcome coefficients + 3 selection coefficients + sigma + rho = 7.
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, start = rep(0.1, 6L)),
               "'start' must have length 7")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, start = rep(0.1, 8L)),
               "'start' must have length 7")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid,
                          start = c(rep(0.1, 5L), NA, 0)),
               "start")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid, start = "auto"),
               "start")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid,
                          start = c(rep(0.1, 5L), 0, 0)),
               "'sigma'.*strictly positive")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid,
                          start = c(rep(0.1, 5L), -1, 0)),
               "'sigma'.*strictly positive")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, valid,
                          start = c(rep(0.1, 5L), 1, 1)),
               "'rho'.*between -1 and 1")
})

test_that("the selection indicator must be a complete 0/1 variable", {
  valid <- make_valid_data()

  wrong_values <- valid
  wrong_values$s[1L] <- 2L
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, wrong_values),
               "only the values 0 and 1.*2")

  fractional <- valid
  fractional$s[1L] <- 0.5
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, fractional),
               "only the values 0 and 1")

  with_na <- valid
  with_na$s[1L] <- NA_integer_
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, with_na),
               "selection indicator.*missing")

  as_factor <- valid
  as_factor$s <- factor(ifelse(valid$s == 1L, "yes", "no"))
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, as_factor),
               "numeric \\(0/1\\) or logical")

  as_character <- valid
  as_character$s <- as.character(valid$s)
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, as_character),
               "numeric \\(0/1\\) or logical")
})

test_that("a logical selection indicator is accepted", {
  valid <- make_valid_data(n = 200L)
  valid$s <- valid$s == 1L
  fit <- poisselect(y ~ x1, s ~ x1 + z1, valid, K = 5L)
  expect_s3_class(fit, "poisselect")
  expect_equal(fit$n_selected, sum(valid$s))
})

test_that("both selection groups must be present", {
  valid <- make_valid_data()

  all_selected <- valid
  all_selected$s <- 1L
  all_selected$y[is.na(all_selected$y)] <- 1L
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, all_selected),
               "All units are selected")

  none_selected <- valid
  none_selected$s <- 0L
  none_selected$y <- NA_integer_
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, none_selected),
               "No unit is selected")
})

test_that("the outcome must be a non-negative integer count", {
  valid <- make_valid_data()
  selected_row <- which(valid$s == 1L)[1L]

  negative <- valid
  negative$y[selected_row] <- -1L
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, negative),
               "non-negative integer counts")

  fractional <- valid
  fractional$y[selected_row] <- 1.5
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, fractional),
               "non-negative integer counts")

  infinite <- valid
  infinite$y[selected_row] <- Inf
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, infinite),
               "non-negative integer counts")

  as_character <- valid
  as_character$y <- as.character(valid$y)
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, as_character),
               "numeric count variable")

  missing_outcome <- valid
  missing_outcome$y[selected_row] <- NA_integer_
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, missing_outcome),
               "outcome is missing for 1 selected unit")

  transformed <- valid
  expect_error(poisselect(log(y + 1) ~ x1, s ~ x1 + z1, transformed),
               "non-negative integer counts")

  all_zero <- valid
  all_zero$y[!is.na(all_zero$y)] <- 0L
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, all_zero),
               "outcome is 0 for every selected unit")
})

test_that("an outcome observed for non-selected units is only reported", {
  valid <- make_valid_data(n = 200L)
  valid$y[is.na(valid$y)] <- 0L
  expect_message(
    fit <- poisselect(y ~ x1, s ~ x1 + z1, valid, K = 5L),
    "observed for .* non-selected unit"
  )
  expect_s3_class(fit, "poisselect")
})

test_that("covariates must be complete where the likelihood uses them", {
  valid <- make_valid_data()
  selected_row <- which(valid$s == 1L)[1L]
  unselected_row <- which(valid$s == 0L)[1L]

  with_na <- valid
  with_na$x2[selected_row] <- NA_real_
  expect_error(poisselect(y ~ x1 + x2, s ~ x1 + z1, with_na),
               "'outcome' equation must be complete for the selected.*'x2'")

  with_inf <- valid
  with_inf$x2[selected_row] <- Inf
  expect_error(poisselect(y ~ x1 + x2, s ~ x1 + z1, with_inf), "'x2'")

  # The selection covariates enter for every unit, including the non-selected
  # ones whose outcome is never used.
  unselected_na <- valid
  unselected_na$z1[unselected_row] <- NA_real_
  expect_error(poisselect(y ~ x1, s ~ x1 + z1, unselected_na),
               "'selection' equation must be complete for all units.*'z1'")
})

test_that("a missing outcome covariate of a non-selected unit is allowed", {
  valid <- make_valid_data(n = 200L)
  valid$x2[which(valid$s == 0L)[1:3]] <- NA_real_
  fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, valid, K = 5L)
  expect_s3_class(fit, "poisselect")
  expect_true(all(is.finite(coef(fit))))
})

test_that("collinear and empty design matrices are rejected by name", {
  valid <- make_valid_data()
  valid$x1_copy <- 2 * valid$x1
  expect_error(poisselect(y ~ x1 + x1_copy, s ~ x1 + z1, valid),
               "'outcome' equation is rank deficient.*'x1_copy'")
  expect_error(poisselect(y ~ x1, s ~ x1 + z1 + x1_copy, valid),
               "'selection' equation is rank deficient.*'x1_copy'")
  # A constant covariate is collinear with the intercept.
  valid$constant <- 1
  expect_error(poisselect(y ~ x1 + constant, s ~ x1 + z1, valid),
               "'constant'")
  expect_error(poisselect(y ~ 0, s ~ x1 + z1, valid),
               "'outcome' equation has no covariates")
})

test_that("the sample must be larger than the number of parameters", {
  small <- make_valid_data(n = 6L)
  expect_error(poisselect(y ~ x1 + x2, s ~ x1 + z1, small),
               "8 parameters but 'data' provides only 6")
  valid <- make_valid_data()
  valid$s <- 0L
  valid$s[which(!is.na(valid$y))[1:3]] <- 1L
  valid$y[valid$s == 0L] <- NA_integer_
  expect_error(poisselect(y ~ x1 + x2, s ~ x1 + z1, valid),
               "observed for only 3 unit")
})

test_that("a missing exclusion restriction produces a warning", {
  # Such a fit is weakly identified and may raise further warnings about the
  # boundary or the Hessian, so only the presence of this one is checked.
  valid <- make_valid_data(n = 200L)
  warnings <- capture_warnings(poisselect(y ~ x1 + x2, s ~ x1 + x2, valid,
                                          K = 5L))
  expect_true(any(grepl("No exclusion restriction", warnings)))
  # The intercept alone does not count as an exclusion restriction.
  warnings <- capture_warnings(poisselect(y ~ 0 + x1 + x2, s ~ x1 + x2, valid,
                                          K = 5L))
  expect_true(any(grepl("No exclusion restriction", warnings)))
  # Extra covariates in the outcome equation only are fine.
  expect_silent(poisselect(y ~ x1 + x2, s ~ x1 + z1, valid, K = 5L))
})

test_that("perfect separation in the selection equation is reported", {
  valid <- make_valid_data(n = 300L, seed = 5L)
  valid$s <- as.integer(valid$z1 > 0)
  valid$y <- ifelse(valid$s == 1L, ifelse(is.na(valid$y), 2L, valid$y), NA)
  warnings <- capture_warnings(poisselect(y ~ x1 + x2, s ~ x1 + z1, valid,
                                          K = 5L))
  expect_true(any(grepl("separate the selected from the non-selected units",
                        warnings)))
})

test_that("factors, interactions and intercept-only equations work", {
  valid <- make_valid_data(n = 300L)
  valid$group <- factor(rep(c("a", "b", "c"), length.out = nrow(valid)))
  fit <- poisselect(y ~ x1 + group, s ~ x1 + z1 + group, valid, K = 5L)
  expect_equal(names(coef(fit, which = "outcome")),
               c("(Intercept)", "x1", "groupb", "groupc"))
  expect_true("groupc" %in% names(coef(fit, which = "selection")))
  fit <- poisselect(y ~ x1 * x2, s ~ x1 + z1, valid, K = 5L)
  expect_true("x1:x2" %in% names(coef(fit, which = "outcome")))
  fit <- poisselect(y ~ 1, s ~ x1 + z1, valid, K = 5L)
  expect_equal(names(coef(fit, which = "outcome")), "(Intercept)")
})

test_that("the checks of simulate_poisselect() catch invalid arguments", {
  expect_error(simulate_poisselect(n = 1L), "n")
  expect_error(simulate_poisselect(n = 10.5), "n")
  expect_error(simulate_poisselect(n = "100"), "n")
  expect_error(simulate_poisselect(n = 100L, beta = c(1, 2)), "beta")
  expect_error(simulate_poisselect(n = 100L, gamma = c(1, 2, 3, 4)), "gamma")
  expect_error(simulate_poisselect(n = 100L, beta = c(1, NA, 3)), "beta")
  expect_error(simulate_poisselect(n = 100L, sigma = 0), "strictly positive")
  expect_error(simulate_poisselect(n = 100L, sigma = -1), "strictly positive")
  expect_error(simulate_poisselect(n = 100L, rho = 1), "between -1 and 1")
  expect_error(simulate_poisselect(n = 100L, rho = -1.5), "between -1 and 1")
  expect_error(simulate_poisselect(n = 100L, rho = NA_real_), "rho")
})

test_that("predict() rejects invalid arguments", {
  expect_error(predict(reference_fit, type = "nonsense"), "should be one of")
  expect_error(predict(reference_fit, newdata = "not a data frame"),
               "newdata")
  expect_error(predict(reference_fit, newdata = c(x1 = 1, x2 = 2)), "newdata")
  expect_error(predict(reference_fit, newdata = data.frame(x1 = numeric(0L))),
               "newdata")
  expect_error(predict(reference_fit, newdata = data.frame(x1 = 1)),
               "'newdata' is missing.*'x2'")
  expect_error(
    predict(reference_fit, newdata = data.frame(x1 = 1), type = "pselect"),
    "'newdata' is missing.*'z1'"
  )
  expect_error(
    predict(reference_fit, newdata = data.frame(x1 = NA_real_, x2 = 1)),
    "'newdata' contains missing values.*'x1'"
  )
})

test_that("plot() and coef() reject invalid arguments", {
  expect_error(plot(reference_fit, which = 3L), "which")
  expect_error(plot(reference_fit, which = 0L), "which")
  expect_error(plot(reference_fit, which = "both"), "which")
  expect_error(plot(reference_fit, which = integer(0L)), "which")
  expect_error(plot(reference_fit, max_count = -1L), "max_count")
  expect_error(plot(reference_fit, max_count = 2.5), "max_count")
  expect_error(plot(reference_fit, n_grid = 2L), "n_grid")
  expect_error(plot(reference_fit, n_grid = NA), "n_grid")
  expect_error(coef(reference_fit, which = "nonsense"), "should be one of")
})
