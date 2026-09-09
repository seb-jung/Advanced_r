test_that("simulate_poisselect() produces data matching the model", {
  set.seed(91L)
  simulated <- simulate_poisselect(n = 5000L, beta = true_beta,
                                   gamma = true_gamma, sigma = true_sigma,
                                   rho = true_rho)
  expect_s3_class(simulated, "data.frame")
  expect_equal(nrow(simulated), 5000L)
  expect_equal(names(simulated), c("y", "s", "x1", "x2", "z1", "y_complete"))
  expect_true(all(simulated$s %in% c(0L, 1L)))
  # y is NA exactly where s = 0, and equals y_complete where s = 1.
  expect_identical(is.na(simulated$y), simulated$s == 0L)
  expect_equal(simulated$y[simulated$s == 1L],
               simulated$y_complete[simulated$s == 1L])
  expect_true(all(simulated$y_complete >= 0L))
  expect_true(all(simulated$y_complete == round(simulated$y_complete)))
  # Share of selected units: s = 1 iff gamma'z + u > 0. Since x1, z1 and u
  # are independent N(0, 1), gamma'z + u is normal with mean gamma_0 = 0.3
  # and variance 1 + 0.5^2 + 0.7^2, so P(s = 1) = Phi(0.3 / sqrt(...)).
  expect_equal(mean(simulated$s), pnorm(0.3 / sqrt(1 + 0.5^2 + 0.7^2)),
               tolerance = 0.05)
  # With rho > 0 the selected units should have larger outcomes on average,
  # that is the selection bias we want to see.
  expect_gt(mean(simulated$y, na.rm = TRUE), mean(simulated$y_complete))
})

test_that("simulate_poisselect() is reproducible under a fixed seed", {
  set.seed(5L)
  first <- simulate_poisselect(n = 100L)
  set.seed(5L)
  second <- simulate_poisselect(n = 100L)
  expect_identical(first, second)
})

test_that("the bundled data set is what its documentation says", {
  expect_equal(dim(sim_selection), c(800L, 6L))
  expect_equal(names(sim_selection), c("y", "s", "x1", "x2", "z1",
                                       "y_complete"))
  expect_identical(is.na(sim_selection$y), sim_selection$s == 0L)
})
