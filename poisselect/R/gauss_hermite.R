#' Gauss-Hermite Nodes and Weights
#'
#' Computes the nodes t_k and weights w_k of the K-point Gauss-Hermite rule
#' for the weight function exp(-t^2). The assignment says we have to compute
#' them ourselves (no quadrature package), so we use the Golub-Welsch method:
#'
#' The monic Hermite polynomials follow the three-term recurrence
#' pi_{k+1}(t) = t * pi_k(t) - (k/2) * pi_{k-1}(t). Writing that recurrence
#' as a matrix gives a symmetric tridiagonal "Jacobi matrix" with zeros on
#' the diagonal and sqrt(k/2) next to it. The eigenvalues of this matrix are
#' exactly the roots of pi_K, i.e. the nodes. The weights are the squared
#' first entries of the normalised eigenvectors, times the total mass of the
#' weight function, which is the integral of exp(-t^2), i.e. sqrt(pi).
#'
#' @param n_nodes Number of quadrature nodes, an integer of at least 2.
#'
#' @return A list with the numeric vectors `nodes` (increasing) and `weights`.
#' @noRd
build_gauss_hermite <- function(n_nodes) {
  jacobi <- matrix(0, n_nodes, n_nodes)
  # Only fill the band above the diagonal, then add the transpose to get the
  # band below as well. Saves us fiddling with two index matrices.
  above <- cbind(seq_len(n_nodes - 1L), seq_len(n_nodes - 1L) + 1L)
  jacobi[above] <- sqrt(seq_len(n_nodes - 1L) / 2)
  jacobi <- jacobi + t(jacobi)
  decomposition <- eigen(jacobi, symmetric = TRUE)
  # eigen() returns the eigenvalues from largest to smallest. We want the
  # nodes in increasing order, so reorder values and vectors together.
  ordering <- order(decomposition$values)
  list(
    nodes = decomposition$values[ordering],
    weights = sqrt(pi) * decomposition$vectors[1L, ordering]^2
  )
}
