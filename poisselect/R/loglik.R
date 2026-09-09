# The approximated log-likelihood from the assignment:
#
#   l(theta) = sum_i { s_i ln[ 1/sqrt(pi) sum_k w_k p(y_i | mu_ik) Phi(eta_ik) ]
#                      + (1 - s_i) ln[1 - Phi(z_i'gamma)] }
#
# with mu_ik = exp(x_i'beta + sqrt(2) sigma t_k)
# and  eta_ik = (z_i'gamma + sqrt(2) rho t_k) / sqrt(1 - rho^2).
#
# Where do the sqrt(2) and the 1/sqrt(pi) come from? The integral in the
# likelihood runs over the standard normal density, but Gauss-Hermite is made
# for the weight function exp(-t^2). Substituting v = sqrt(2) * t turns
# exp(-v^2 / 2) into exp(-t^2); what is left over from dv / sqrt(2 pi) is
# exactly dt / sqrt(pi). That is why every node enters as sqrt(2) * t_k and
# the whole sum is divided by sqrt(pi).
#
# The second important point: everything is done on the log scale. For a
# count like y = 500 the Poisson probability is around 1e-300 and the product
# with the weights would underflow to 0 (and log(0) = -Inf). So we never
# build the inner sum on the probability scale, see log_sum_exp_rows().

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
  # Non-selected units only contribute ln(1 - Phi(z'gamma)). We must NOT
  # write log(1 - pnorm(...)) here: for z'gamma = 10, pnorm() is already
  # exactly 1 in double precision, so 1 - pnorm() = 0 and the log is -Inf.
  # pnorm() can compute the upper tail and its log directly without that
  # cancellation problem.
  loglik_unselected <- pnorm(
    drop(model$z_unselected %*% parameters$gamma),
    lower.tail = FALSE,
    log.p = TRUE
  )
  sum(log_sum_exp_rows(log_terms)) + sum(loglik_unselected)
}

#' Node-Specific Quantities of the Selected Units
#'
#' Builds the n_selected x K matrices that show up in the likelihood. The
#' gradient and plot 1 need exactly the same matrices, so they are computed
#' in one place (DRY) and only combined differently afterwards.
#'
#' Why eta looks like that: given the standardised outcome error v, the
#' selection error u is conditionally normal with mean rho * v and variance
#' 1 - rho^2. So Pr(u > -z'gamma | v) = Phi((z'gamma + rho v) /
#' sqrt(1 - rho^2)), and at node k we plug in v = sqrt(2) * t_k.
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
  # outer() builds the whole n_selected x K matrix in one go: row i is the
  # linear predictor of unit i plus each of the K scaled nodes. No loop over
  # observations needed (vectorisation, see the performance lecture).
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
  # dpois() with a vector y and a matrix of rates: R recycles y down the
  # columns (column-major storage), so row i really gets y_i in every column.
  # dpois() returns a plain vector though, so we put the dim back on.
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
#' The log-sum-exp trick. We want log(sum_k exp(a_ik)), but the a_ik can be
#' very negative (like -700), so exp() would underflow to 0. Instead we take
#' the row maximum A_i out first: log(sum_k exp(a_ik)) = A_i +
#' log(sum_k exp(a_ik - A_i)). Now the largest term is exp(0) = 1 and nothing
#' can underflow. Mathematically identical, numerically safe.
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
#' Analytic gradient of the log-likelihood. optim() can also work with
#' numerical differences, but we found that BFGS then stops too early when
#' the counts get large (up to 47 log-likelihood units below the real
#' optimum). With the analytic gradient that gap is gone and the fit is
#' about 3-4 times faster.
#'
#' How it is derived: for a selected unit l_i = ln sum_k exp(a_ik). The
#' derivative of a log-sum-exp is the weighted average of the derivatives of
#' the single terms, with weights p_ik = exp(a_ik - l_i) (they sum to 1 per
#' row). The single terms are
#'   - a Poisson log-density: derivative (y_i - mu_ik) times d(ln mu_ik),
#'   - ln Phi(eta_ik): derivative phi(eta)/Phi(eta) (the inverse Mills
#'     ratio) times d(eta_ik).
#' Non-selected units only contribute the usual probit gradient of
#' ln(1 - Phi(z'gamma)).
#'
#' Because optim() works on log(sigma) and atanh(rho), the last two entries
#' need the chain rule: d sigma / d log(sigma) = sigma and
#' d rho / d atanh(rho) = 1 - rho^2. Checked against central differences in
#' test-loglik.R.
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
  # Inverse Mills ratio phi(eta) / Phi(eta). Computed as exp(log phi - log
  # Phi), because for very negative eta both would be 0 and 0 / 0 = NaN.
  mills <- exp(dnorm(pieces$eta, log = TRUE) - pieces$log_selection)
  scaled_nodes <- matrix(sqrt(2) * quadrature$nodes, nrow = nrow(log_terms),
                         ncol = ncol(log_terms), byrow = TRUE)
  # Derivative of eta with respect to rho (quotient rule on the formula in
  # compute_node_pieces): (sqrt(2) t_k + rho z_i'gamma) / (1 - rho^2)^(3/2).
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
#' optim() minimises, so we hand it the negative log-likelihood. If the value
#' is not finite (happens when the line search jumps to a crazy region where
#' every node underflows) we return a big finite number instead of Inf/NaN.
#' optim() would abort on NaN; with the penalty it just steps back.
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
