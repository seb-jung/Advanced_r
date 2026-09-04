# The approximated log-likelihood of the assignment,
#
#   l(theta) = sum_i { s_i ln[ 1/sqrt(pi) sum_k w_k p(y_i | mu_ik) Phi(eta_ik) ]
#                      + (1 - s_i) ln[1 - Phi(z_i'gamma)] },
#
# with mu_ik = exp(x_i'beta + sqrt(2) sigma t_k) and
# eta_ik = (z_i'gamma + sqrt(2) rho t_k) / sqrt(1 - rho^2). The factor sqrt(2)
# comes from the substitution v = sqrt(2) t that turns the standard normal
# density into the Gauss-Hermite weight function exp(-t^2), and 1/sqrt(pi) is
# the remaining normalising constant of that substitution.
#
# Everything is evaluated on the log scale: the inner sum is never formed on
# the probability scale, because p(y | mu) underflows for large counts.

#' Approximated Log-Likelihood
#'
#' @param theta Numeric vector of parameters on the unconstrained scale, see
#'   [split_parameters()].
#' @param model A model list as built by [build_model_data()].
#' @param quadrature A list with the components `nodes` and `weights`.
#'
#' @return The scalar value of the approximated log-likelihood.
#' @noRd
compute_loglik <- function(theta, model, quadrature) {
  parameters <- split_parameters(theta, model$n_beta, model$n_gamma)
  pieces <- compute_node_pieces(parameters, model, quadrature)
  log_terms <- compute_log_terms(model$y_selected, pieces)
  # Non-selected units contribute ln(1 - Phi(z'gamma)); lower.tail = FALSE
  # evaluates that upper tail directly, whereas log(1 - pnorm(.)) would lose
  # all precision for large z'gamma.
  loglik_unselected <- pnorm(
    drop(model$z_unselected %*% parameters$gamma),
    lower.tail = FALSE,
    log.p = TRUE
  )
  sum(log_sum_exp_rows(log_terms)) + sum(loglik_unselected)
}

#' Node-Specific Quantities of the Selected Units
#'
#' Shared by the log-likelihood, its gradient and the model-implied count
#' distribution of [plot.poisselect()], which combine the same pieces in
#' different ways.
#'
#' @param parameters A parameter list as returned by [split_parameters()].
#' @param model A model list as built by [build_model_data()].
#' @param quadrature A list with the components `nodes` and `weights`.
#'
#' @return A list with the `n_selected` by `K` matrices `log_rate`
#'   (\eqn{\ln \mu_{ik}}), `eta` and `log_selection` (\eqn{\ln\Phi(\eta_{ik})})
#'   and the length-`K` vector `log_weight` (\eqn{\ln w_k - \ln\pi / 2}).
#' @noRd
compute_node_pieces <- function(parameters, model, quadrature) {
  scaled_nodes <- sqrt(2) * quadrature$nodes
  # outer() adds the K scaled nodes to every unit's linear predictor and gives
  # the n_selected x K matrices in one step.
  log_rate <- outer(drop(model$x_selected %*% parameters$beta),
                    parameters$sigma * scaled_nodes, "+")
  eta <- outer(drop(model$z_selected %*% parameters$gamma),
               parameters$rho * scaled_nodes, "+") / sqrt(1 - parameters$rho^2)
  list(
    log_rate = log_rate,
    eta = eta,
    log_selection = pnorm(eta, log.p = TRUE),
    log_weight = log(quadrature$weights) - 0.5 * log(pi)
  )
}

#' Logarithmic Summands a_ik of the Inner Sum
#'
#' @param y Numeric vector of counts, one per selected unit, or a single count.
#' @param pieces A list as returned by [compute_node_pieces()].
#'
#' @return The `n_selected` by `K` matrix of
#'   \eqn{a_{ik} = \ln w_k - \ln\pi/2 + \ln p(y_i \mid \mu_{ik}) +
#'   \ln\Phi(\eta_{ik})}.
#' @noRd
compute_log_terms <- function(y, pieces) {
  # y is recycled down the columns of the matrix, so every row keeps its own
  # count; the dimension is restored explicitly because dpois() drops it.
  log_density <- dpois(y, exp(pieces$log_rate), log = TRUE)
  dim(log_density) <- dim(pieces$log_rate)
  add_log_weight(log_density + pieces$log_selection, pieces$log_weight)
}

#' Add the Logarithmic Quadrature Weight to Every Column
#'
#' @param a Numeric matrix with one column per quadrature node.
#' @param log_weight Numeric vector of logarithmic weights, one per node.
#'
#' @return Numeric matrix of the same dimension as `a`.
#' @noRd
add_log_weight <- function(a, log_weight) {
  a + matrix(log_weight, nrow = nrow(a), ncol = ncol(a), byrow = TRUE)
}

#' Row-Wise Log-Sum-Exp
#'
#' Subtracting the row maximum before exponentiating keeps the largest summand
#' at exp(0) = 1, so the tiny Poisson probabilities of large counts cannot
#' underflow to zero.
#'
#' @param a Numeric matrix of logarithmic summands, one row per unit.
#'
#' @return Numeric vector with `log(rowSums(exp(a)))`, one entry per row.
#' @noRd
log_sum_exp_rows <- function(a) {
  row_maximum <- apply(a, 1L, max)
  row_maximum + log(rowSums(exp(a - row_maximum)))
}

#' Gradient of the Approximated Log-Likelihood
#'
#' For a selected unit the contribution is \eqn{\ell_i = \ln\sum_k e^{a_{ik}}},
#' whose derivative is the average of the node derivatives
#' \eqn{\partial a_{ik}/\partial\theta} weighted by
#' \eqn{p_{ik} = e^{a_{ik} - \ell_i}}. The node derivatives are those of a
#' Poisson log-density, \eqn{(y_i - \mu_{ik})} times the derivative of
#' \eqn{\ln\mu_{ik}}, and of \eqn{\ln\Phi(\eta_{ik})}, which is the inverse
#' Mills ratio \eqn{\phi(\eta)/\Phi(\eta)} times the derivative of
#' \eqn{\eta_{ik}}. The non-selected units only add the probit gradient of
#' \eqn{\ln[1 - \Phi(z_i'\gamma)]}. The chain rule for the unconstrained scale
#' multiplies the `sigma` entry by \eqn{d\sigma/d\ln\sigma = \sigma} and the
#' `rho` entry by \eqn{d\rho/d\,\mathrm{atanh}\rho = 1 - \rho^2}.
#'
#' @param theta Numeric vector of parameters on the unconstrained scale.
#' @param model A model list as built by [build_model_data()].
#' @param quadrature A list with the components `nodes` and `weights`.
#'
#' @return Numeric vector of the same length as `theta`.
#' @noRd
compute_loglik_gradient <- function(theta, model, quadrature) {
  parameters <- split_parameters(theta, model$n_beta, model$n_gamma)
  pieces <- compute_node_pieces(parameters, model, quadrature)
  log_terms <- compute_log_terms(model$y_selected, pieces)
  node_weight <- exp(log_terms - log_sum_exp_rows(log_terms))
  residual <- model$y_selected - exp(pieces$log_rate)
  # phi / Phi evaluated on the log scale, so extreme eta cannot give 0 / 0.
  mills <- exp(dnorm(pieces$eta, log = TRUE) - pieces$log_selection)
  scaled_nodes <- matrix(sqrt(2) * quadrature$nodes, nrow = nrow(log_terms),
                         ncol = ncol(log_terms), byrow = TRUE)
  # d eta_ik / d rho = (sqrt(2) t_k + rho z_i'gamma) / (1 - rho^2)^(3/2).
  linear_selected <- drop(model$z_selected %*% parameters$gamma)
  d_eta_d_rho <- (scaled_nodes + parameters$rho * linear_selected) /
    (1 - parameters$rho^2)^1.5
  linear_unselected <- drop(model$z_unselected %*% parameters$gamma)
  mills_unselected <- exp(
    dnorm(linear_unselected, log = TRUE) -
      pnorm(linear_unselected, lower.tail = FALSE, log.p = TRUE)
  )
  d_beta <- crossprod(model$x_selected, rowSums(node_weight * residual))
  d_gamma <- crossprod(model$z_selected, rowSums(node_weight * mills)) /
    sqrt(1 - parameters$rho^2) -
    crossprod(model$z_unselected, mills_unselected)
  d_sigma <- sum(node_weight * residual * scaled_nodes)
  d_rho <- sum(node_weight * mills * d_eta_d_rho)
  c(d_beta, d_gamma, d_sigma * parameters$sigma, d_rho * (1 - parameters$rho^2))
}

#' Negative Log-Likelihood and Gradient for the Optimiser
#'
#' [stats::optim()] minimises, so the sign is flipped. A non-finite value is
#' replaced by a large finite penalty: it only occurs when the optimiser probes
#' a region in which every quadrature node underflows, and the penalty pushes
#' the line search back instead of aborting the whole fit.
#'
#' @inheritParams compute_loglik
#'
#' @return A scalar, respectively a numeric vector.
#' @noRd
compute_negative_loglik <- function(theta, model, quadrature) {
  value <- -compute_loglik(theta, model, quadrature)
  if (!is.finite(value)) {
    return(1e10)
  }
  value
}

#' @rdname compute_negative_loglik
#' @noRd
compute_negative_gradient <- function(theta, model, quadrature) {
  -compute_loglik_gradient(theta, model, quadrature)
}
