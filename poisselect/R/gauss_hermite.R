#' Gauss-Hermite Nodes and Weights
#'
#' Computes the nodes \eqn{t_k} and weights \eqn{w_k} of the `n_nodes`-point
#' Gauss-Hermite rule for the weight function \eqn{e^{-t^2}} with the
#' Golub-Welsch method. The monic Hermite polynomials satisfy the recurrence
#' \eqn{\pi_{k+1}(t) = t\,\pi_k(t) - (k/2)\,\pi_{k-1}(t)}, so their Jacobi
#' matrix is symmetric tridiagonal with zeros on the diagonal and
#' \eqn{\sqrt{k/2}} next to it. Its eigenvalues are the nodes, and the weights
#' are the squared first components of the normalised eigenvectors times
#' \eqn{\mu_0 = \int e^{-t^2}\,dt = \sqrt{\pi}}.
#'
#' @param n_nodes Number of quadrature nodes, an integer of at least 2.
#'
#' @return A list with the numeric vectors `nodes` (increasing) and `weights`.
#' @noRd
build_gauss_hermite <- function(n_nodes) {
  jacobi <- matrix(0, n_nodes, n_nodes)
  # Fill the band above the diagonal, then mirror it by adding the transpose.
  above <- cbind(seq_len(n_nodes - 1L), seq_len(n_nodes - 1L) + 1L)
  jacobi[above] <- sqrt(seq_len(n_nodes - 1L) / 2)
  jacobi <- jacobi + t(jacobi)
  decomposition <- eigen(jacobi, symmetric = TRUE)
  # eigen() sorts the eigenvalues in decreasing order; the rule is reported
  # with increasing nodes.
  ordering <- order(decomposition$values)
  list(
    nodes = decomposition$values[ordering],
    weights = sqrt(pi) * decomposition$vectors[1L, ordering]^2
  )
}
