#' Translate the Two Model Formulas into Responses and Design Matrices
#'
#' Turns the two formulas plus the data.frame into everything the likelihood
#' needs. Both equations go through the same helper (DRY), so model frames,
#' design matrices, `terms` and factor levels are all built the same way.
#'
#' We also split the rows into selected / non-selected right here, once.
#' optim() calls the likelihood hundreds of times and it would be wasteful
#' to redo the subsetting in every single call.
#'
#' @param outcome Two-sided formula of the outcome equation.
#' @param selection Two-sided formula of the selection equation.
#' @param data A `data.frame` containing all variables of both formulas.
#'
#' @return A list with the responses `y` and `s`, the logical vector
#'   `selected`, the full design matrices `x` and `z`, their selected parts
#'   `x_selected`, `z_selected`, `y_selected`, the non-selected part
#'   `z_unselected`, the lists `terms` and `xlevels` (one entry per equation)
#'   and the counts `n`, `n_selected`, `n_beta` and `n_gamma`.
#' @noRd
build_model_data <- function(outcome, selection, data) {
  parts <- lapply(
    list(outcome = outcome, selection = selection),
    build_equation_data,
    data = data
  )
  s <- as_selection_indicator(parts$selection$response)
  selected <- !is.na(s) & s == 1
  x <- parts$outcome$design
  z <- parts$selection$design
  list(
    y = parts$outcome$response,
    s = s,
    selected = selected,
    x = x,
    z = z,
    x_selected = x[selected, , drop = FALSE],
    z_selected = z[selected, , drop = FALSE],
    z_unselected = z[!selected, , drop = FALSE],
    y_selected = parts$outcome$response[selected],
    terms = lapply(parts, `[[`, "terms"),
    xlevels = lapply(parts, `[[`, "xlevels"),
    n = length(s),
    n_selected = sum(selected),
    n_beta = ncol(x),
    n_gamma = ncol(z)
  )
}

#' Model Frame, Design Matrix and Terms of a Single Equation
#'
#' Important detail: `na.action = na.pass`. The default na.omit would throw
#' away every row with an NA, but y is NA for all non-selected units by
#' design. So we keep the NAs here and let [check_model_data()] decide
#' afterwards which NAs are fine and which ones are errors.
#'
#' @param formula A two-sided formula.
#' @param data A `data.frame`.
#'
#' @return A list with `response`, `design`, `terms` and `xlevels`.
#' @noRd
build_equation_data <- function(formula, data) {
  model_terms <- terms(formula, data = data)
  frame <- model.frame(model_terms, data = data, na.action = na.pass)
  list(
    response = model.response(frame),
    design = model.matrix(model_terms, frame),
    terms = model_terms,
    xlevels = .getXlevels(model_terms, frame)
  )
}

#' Coerce the Selection Response to a Numeric 0/1 Indicator
#'
#' TRUE/FALSE is accepted and turned into 1/0. Factors and character vectors
#' are rejected: we could map the levels to 0/1, but which level would be
#' "selected"? It would depend on the alphabetical order of the labels, and
#' a silent wrong guess flips the whole selection equation. So the user has
#' to convert explicitly.
#'
#' @param response The response of the selection equation.
#'
#' @return A numeric vector.
#' @noRd
as_selection_indicator <- function(response) {
  if (is.logical(response)) {
    return(as.numeric(response))
  }
  if (!is.numeric(response)) {
    stop("The response of 'selection' must be numeric (0/1) or logical, but ",
         "is of class ", sQuote(class(response)[1L]), ". Please convert the ",
         "indicator explicitly, for example with 'as.integer(s == \"yes\")'.",
         call. = FALSE)
  }
  as.numeric(response)
}
