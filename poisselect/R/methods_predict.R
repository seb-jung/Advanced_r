#' Predictions from a Poisson Selection Model
#'
#' @details
#' The three types answer different questions.
#'
#' * `"link"` returns the linear predictor \eqn{x'\hat\beta} of the outcome
#'   equation.
#' * `"response"` returns the *unconditional* population mean
#'   \eqn{\hat E[Y \mid x] = \exp(x'\hat\beta + \hat\sigma^2/2)}, that is the
#'   mean averaged over the unobserved outcome heterogeneity
#'   \eqn{\varepsilon}. This is deliberately not the mean of the selected
#'   subpopulation \eqn{E[Y \mid x, s = 1]}, which would still carry the
#'   selection bias. The factor \eqn{e^{\hat\sigma^2/2}} is the mean of the
#'   log-normal multiplicative error, so `"response"` exceeds `exp("link")`.
#' * `"pselect"` returns the selection probability \eqn{\Phi(z'\hat\gamma)}.
#'
#' @param object An object of class `"poisselect"`.
#' @param newdata Optional `data.frame` with the covariates at which to
#'   predict. Defaults to `NULL`, which reuses the data the model was fitted
#'   on. `"link"` and `"response"` need the covariates of the outcome equation,
#'   `"pselect"` those of the selection equation.
#' @param type Type of prediction, one of `"response"` (the default), `"link"`
#'   or `"pselect"`.
#' @param ... Currently ignored, present for consistency with the generic.
#'
#' @return A numeric vector with one entry per row of `newdata`, or per row of
#'   the original data if `newdata` is `NULL`.
#'
#' @examples
#' fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)
#'
#' head(predict(fit))
#' head(predict(fit, type = "link"))
#' head(predict(fit, type = "pselect"))
#'
#' # Predictions for new covariate values.
#' grid <- data.frame(x1 = c(-1, 0, 1), x2 = 0, z1 = 0)
#' predict(fit, newdata = grid)
#' predict(fit, newdata = grid, type = "pselect")
#'
#' @export
predict.poisselect <- function(object, newdata = NULL,
                               type = c("response", "link", "pselect"), ...) {
  type <- match.arg(type)
  assert_data_frame(newdata, min.rows = 1L, null.ok = TRUE,
                    .var.name = "newdata")
  equation <- switch(type, pselect = "selection", "outcome")
  design <- build_prediction_matrix(object, newdata, equation)
  linear_predictor <- drop(design %*% object$coefficients[[equation]])
  switch(
    type,
    link = linear_predictor,
    response = exp(linear_predictor + object$sigma^2 / 2),
    pselect = pnorm(linear_predictor)
  )
}

#' Design Matrix for a Prediction
#'
#' Without newdata we simply return the design matrix stored in the fit.
#' With newdata we rebuild it from the stored `terms` and factor levels
#' (xlevels), the same way predict.lm() does it, so that the dummy coding of
#' factors is identical to the fit even if newdata has fewer levels.
#'
#' @param object An object of class `"poisselect"`.
#' @param newdata A `data.frame`, or `NULL` to reuse the fitted data.
#' @param equation Either `"outcome"` or `"selection"`.
#'
#' @return A numeric design matrix.
#' @noRd
build_prediction_matrix <- function(object, newdata, equation) {
  if (is.null(newdata)) {
    return(object$model[[switch(equation, outcome = "x", selection = "z")]])
  }
  model_terms <- delete.response(object$model$terms[[equation]])
  missing_variables <- setdiff(all.vars(model_terms), names(newdata))
  if (length(missing_variables) > 0L) {
    stop("'newdata' is missing the following variable(s) of the '", equation,
         "' equation: ", toString(sQuote(missing_variables)), ".",
         call. = FALSE)
  }
  frame <- model.frame(model_terms, data = newdata, na.action = na.pass,
                       xlev = object$model$xlevels[[equation]])
  # Check for NA on the model frame, not on the model matrix: model.matrix()
  # would turn an NA column into a weird dummy called e.g. "x1TRUE" and the
  # error message would name that instead of the real variable.
  incomplete <- names(frame)[colSums(is.na(frame)) > 0L]
  if (length(incomplete) > 0L) {
    stop("'newdata' contains missing values in the variable(s) ",
         toString(sQuote(incomplete)), ".", call. = FALSE)
  }
  model.matrix(model_terms, frame)
}
