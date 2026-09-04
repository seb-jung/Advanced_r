# The rule is checked against published values, against the exact moments of
# the weight function and against the way the likelihood uses it.

test_that("the rule reproduces the published values for K = 2, 3 and 5", {
  two_point <- build_gauss_hermite(2L)
  expect_equal(two_point$nodes, c(-1, 1) / sqrt(2))
  expect_equal(two_point$weights, rep(sqrt(pi) / 2, 2L))

  three_point <- build_gauss_hermite(3L)
  expect_equal(three_point$nodes, c(-sqrt(3 / 2), 0, sqrt(3 / 2)))
  expect_equal(three_point$weights,
               c(sqrt(pi) / 6, 2 * sqrt(pi) / 3, sqrt(pi) / 6))

  five_point <- build_gauss_hermite(5L)
  expect_equal(
    five_point$nodes,
    c(-2.020182870456086, -0.958572464613819, 0,
      0.958572464613819, 2.020182870456086),
    tolerance = 1e-12
  )
  expect_equal(
    five_point$weights,
    c(0.019953242059046, 0.393619323152241, 0.945308720482942,
      0.393619323152241, 0.019953242059046),
    tolerance = 1e-12
  )
})

test_that("the weights integrate the weight function exactly", {
  for (n_nodes in c(2L, 5L, 20L, 50L, 200L)) {
    expect_equal(sum(build_gauss_hermite(n_nodes)$weights), sqrt(pi))
  }
})

test_that("a K-point rule is exact for polynomials up to degree 2K - 1", {
  # Odd moments of exp(-t^2) vanish, the even ones are Gamma(m + 1/2). The
  # error is measured relative to the size of the summands, which is the
  # precision floating point arithmetic can deliver for high degrees.
  for (n_nodes in c(2L, 3L, 5L, 10L)) {
    quadrature <- build_gauss_hermite(n_nodes)
    for (degree in 0:(2L * n_nodes - 1L)) {
      analytic <- if (degree %% 2L == 1L) 0 else gamma((degree + 1) / 2)
      terms <- quadrature$weights * quadrature$nodes^degree
      expect_lt(abs(sum(terms) - analytic) / sum(abs(terms)), 1e-12)
    }
    # Degree 2K is the first one the rule does not integrate exactly.
    degree <- 2L * n_nodes
    expect_gt(
      abs(sum(quadrature$weights * quadrature$nodes^degree) -
            gamma((degree + 1) / 2)),
      1e-6
    )
  }
})

test_that("nodes are sorted and symmetric around zero", {
  quadrature <- build_gauss_hermite(20L)
  expect_false(is.unsorted(quadrature$nodes))
  expect_equal(quadrature$nodes, -rev(quadrature$nodes))
  expect_equal(quadrature$weights, rev(quadrature$weights))
  expect_true(all(quadrature$weights > 0))
})

test_that("the rule integrates a standard normal expectation", {
  # With v = sqrt(2) t the rule integrates against the standard normal density,
  # which is exactly how the likelihood uses it: E[exp(v)] = exp(1/2).
  quadrature <- build_gauss_hermite(20L)
  standardised <- quadrature$weights / sqrt(pi)
  expect_equal(sum(standardised * exp(sqrt(2) * quadrature$nodes)), exp(0.5),
               tolerance = 1e-8)
  expect_equal(sum(standardised * (sqrt(2) * quadrature$nodes)^2), 1,
               tolerance = 1e-10)
})
