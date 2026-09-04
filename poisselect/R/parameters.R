# The model requires sigma > 0 and -1 < rho < 1, but optim() works best
# without restrictions. The optimisation therefore runs on the unconstrained
# vector
#
#   theta = [ beta (n_beta), gamma (n_gamma), log(sigma), atanh(rho) ],
#
# whose back-transformations exp() and tanh() map the whole real line into the
# admissible ranges. Fisher's z transformation atanh(rho) is the standard
# choice for a correlation.

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
#' Both equations usually share names, most obviously the intercept, so the
#' equation is used as a prefix.
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
