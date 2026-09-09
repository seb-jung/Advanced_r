# Tests for our Gauss-Hermite implementation. Since we compute the nodes and
# weights ourselves we compare them with published tables (K = 2, 3, 5), with
# the exact moments of exp(-t^2), and with the way the likelihood uses them.

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
  # A K-point Gauss rule integrates every polynomial up to degree 2K - 1
  # exactly. The moments of exp(-t^2) are known: odd ones are 0, even ones
  # are Gamma(m + 1/2). We measure the error relative to the size of the
  # summands, because for high degrees the terms get large and an absolute
  # 1e-12 would be stricter than floating point can deliver.
  for (n_nodes in c(2L, 3L, 5L, 10L)) {
    quadrature <- build_gauss_hermite(n_nodes)
    for (degree in 0:(2L * n_nodes - 1L)) {
      analytic <- if (degree %% 2L == 1L) 0 else gamma((degree + 1) / 2)
      terms <- quadrature$weights * quadrature$nodes^degree
      expect_lt(abs(sum(terms) - analytic) / sum(abs(terms)), 1e-12)
    }
    # And degree 2K should be the first one that is NOT exact anymore, which
    # confirms we really have a K-point rule and not something else.
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
  # Same substitution as in the likelihood: with v = sqrt(2) * t and weights
  # w / sqrt(pi) the rule integrates against the standard normal density.
  # So E[exp(v)] should be exp(1/2) (log-normal mean) and E[v^2] = 1.
  quadrature <- build_gauss_hermite(20L)
  standardised <- quadrature$weights / sqrt(pi)
  expect_equal(sum(standardised * exp(sqrt(2) * quadrature$nodes)), exp(0.5),
               tolerance = 1e-8)
  expect_equal(sum(standardised * (sqrt(2) * quadrature$nodes)^2), 1,
               tolerance = 1e-10)
})
