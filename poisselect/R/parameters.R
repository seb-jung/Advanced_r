# The model needs sigma > 0 and -1 < rho < 1, but BFGS in optim() does not
# know about bounds. Our solution: optimise over an unconstrained vector
#
#   theta = [ beta (n_beta), gamma (n_gamma), log(sigma), atanh(rho) ]
#
# and transform back with exp() and tanh(). Both map the whole real line
# into the allowed range, so the restrictions hold automatically and optim()
# can never produce something like sigma = -0.3 or rho = 1.2. atanh(rho) is
# Fisher's z-transformation, the usual choice for correlations.
#
# The functions below convert between theta and the named parameters, so
# nobody else in the package has to know the layout of theta.

#' Split the Unconstrained Vector into Model Parameters
#'
#' @param theta Numeric vector on the unconstrained scale.
#' @param n_beta Number of outcome coefficients.
#' @param n_gamma Number of selection coefficients.
#'
#' @return A list with `beta`, `gamma`, `sigma` and `rho` on the original scale.
#' @noRd
split_parameters <- function(theta, n_beta, n_gamma) {
  list(
    beta = theta[seq_len(n_beta)],
    gamma = theta[n_beta + seq_len(n_gamma)],
    sigma = exp(theta[n_beta + n_gamma + 1L]),
    rho = tanh(theta[n_beta + n_gamma + 2L])
  )
}

#' Collect Model Parameters into the Unconstrained Vector
#'
#' @param beta,gamma Numeric coefficient vectors.
#' @param sigma Positive scalar.
#' @param rho Scalar in `(-1, 1)`.
#'
#' @return Numeric vector on the unconstrained scale.
#' @noRd
pack_parameters <- function(beta, gamma, sigma, rho) {
  unname(c(beta, gamma, log(sigma), atanh(rho)))
}

#' Names of the Full Parameter Vector
#'
#' Both equations usually have some names in common (at least the
#' intercept, often x1 too), so we prefix them with the equation name.
#' We use "_" and not ":" because in R a ":" in a name looks like an
#' interaction term.
#'
#' @param model A model list as built by [build_model_data()].
#'
#' @return Character vector of length `n_beta + n_gamma + 2`.
#' @noRd
build_parameter_names <- function(model) {
  c(paste0("outcome_", colnames(model$x)),
    paste0("selection_", colnames(model$z)),
    "sigma", "rho")
}
