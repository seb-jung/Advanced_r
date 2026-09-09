# All the input checks for poisselect(). We split them into two stages, like
# the "fail fast, fail appropriately" idea from the lecture:
#
#   1. check_arguments(): looks only at what the user typed (formulas, data,
#      K, start, control). Runs before we touch the data at all, so obvious
#      mistakes fail immediately.
#   2. check_model_data(): runs after the design matrices are built. Some
#      problems (NA pattern, collinear columns, missing exclusion
#      restriction) can only be seen at that point.
#
# Every error message tries to say which argument is wrong, what is wrong
# with it and, where possible, which variable is affected.

#' Check the Arguments of poisselect()
#'
#' @param outcome,selection The two model formulas.
#' @param data The data set.
#' @param n_nodes The number of quadrature nodes.
#' @param start Optional numeric vector of starting values.
#' @param control A list of control parameters for [stats::optim()].
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_arguments <- function(outcome, selection, data, n_nodes, start,
                            control) {
  if (missing(data)) {
    stop("'data' is missing. Please pass a data.frame that contains all ",
         "variables of both formulas.", call. = FALSE)
  }
  check_formula(outcome, "outcome", "y ~ x1 + x2")
  check_formula(selection, "selection", "s ~ x1 + z1")
  assert_data_frame(data, min.rows = 1L, min.cols = 1L, .var.name = "data")
  check_variables_available(list(outcome, selection), data)
  # K = 1 is not allowed: the only node would sit at t = 0, so sigma would
  # drop out of the likelihood completely. The upper limit of 200 is a bit
  # arbitrary, but above ~50 nodes nothing improves anymore and every
  # likelihood evaluation just gets slower.
  assert_int(n_nodes, lower = 2L, upper = 200L, .var.name = "K")
  assert_numeric(start, any.missing = FALSE, finite = TRUE, null.ok = TRUE,
                 .var.name = "start")
  assert_list(control, .var.name = "control")
  if (length(control) > 0L) {
    assert_names(
      names(control),
      subset.of = c("maxit", "reltol", "abstol", "trace", "REPORT", "fnscale",
                    "parscale", "ndeps"),
      .var.name = "names(control)"
    )
  }
  invisible(TRUE)
}

#' Check a Single Model Formula
#'
#' @param formula The object to check.
#' @param name The argument name, used in the error messages.
#' @param example A valid formula shown in the error message.
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_formula <- function(formula, name, example) {
  assert_formula(formula, .var.name = name)
  if (length(formula) != 3L) {
    stop("'", name, "' must be a two-sided formula such as '", example,
         "', but no left-hand side was given.", call. = FALSE)
  }
  # We forbid the '.' shortcut on purpose. Both formulas use the same data
  # set, so 'y ~ .' would silently put s (and z1 etc.) into the outcome
  # equation, which makes no sense.
  if ("." %in% all.vars(formula)) {
    stop("'", name, "' must not use the '.' shortcut, because both equations ",
         "are built from the same data set. Please list the covariates ",
         "explicitly, as in '", example, "'.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check that All Formula Variables Exist in the Data
#'
#' @param formulas A list of formulas.
#' @param data A `data.frame`.
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_variables_available <- function(formulas, data) {
  required <- unique(unlist(lapply(formulas, all.vars)))
  missing_variables <- setdiff(required, names(data))
  if (length(missing_variables) > 0L) {
    stop("The following variable(s) are used in the formulas but are not ",
         "columns of 'data': ", toString(sQuote(missing_variables)), ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Check the User-Supplied Starting Values
#'
#' This check cannot be part of check_arguments(), because we only know how
#' long `start` has to be after the design matrices are built (number of
#' columns of x and z, plus 2 for sigma and rho).
#'
#' @param start The starting values, or `NULL`.
#' @param model A model list as built by [build_model_data()].
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_start_values <- function(start, model) {
  if (is.null(start)) {
    return(invisible(TRUE))
  }
  n_parameter <- model$n_beta + model$n_gamma + 2L
  if (length(start) != n_parameter) {
    stop("'start' must have length ", n_parameter, ": ", model$n_beta,
         " outcome coefficients, ", model$n_gamma, " selection coefficients, ",
         "then sigma and rho, but has length ", length(start), ".",
         call. = FALSE)
  }
  sigma_start <- start[n_parameter - 1L]
  rho_start <- start[n_parameter]
  if (sigma_start <= 0) {
    stop("The starting value for 'sigma' (element ", n_parameter - 1L,
         " of 'start') must be strictly positive, but is ", sigma_start, ".",
         call. = FALSE)
  }
  if (abs(rho_start) >= 1) {
    stop("The starting value for 'rho' (element ", n_parameter, " of 'start') ",
         "must lie strictly between -1 and 1, but is ", rho_start, ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Check the Data Implied by the Two Formulas
#'
#' @param model A model list as built by [build_model_data()].
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_model_data <- function(model) {
  check_selection_indicator(model$s)
  check_outcome_counts(model$y, model$selected)
  check_sample_size(model)
  # Which rows need to be complete depends on the equation: the outcome
  # covariates are only used for the selected units (y is unknown for the
  # others anyway), but the selection covariates are needed for everyone.
  check_design_matrix(model$x_selected, "outcome")
  check_design_matrix(model$z, "selection")
  check_exclusion_restriction(model)
  invisible(TRUE)
}

#' Check the Selection Indicator
#'
#' @param s The selection indicator, already coerced to numeric.
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_selection_indicator <- function(s) {
  if (anyNA(s)) {
    stop("The selection indicator (left-hand side of 'selection') contains ",
         sum(is.na(s)), " missing value(s), but it must be observed for ",
         "every unit.", call. = FALSE)
  }
  if (!all(s %in% c(0, 1))) {
    stop("The selection indicator (left-hand side of 'selection') must ",
         "contain only the values 0 and 1, but also contains ",
         toString(head(sort(unique(s[!s %in% c(0, 1)])), 5L)), ".",
         call. = FALSE)
  }
  if (all(s == 1)) {
    stop("All units are selected. Without non-selected units the selection ",
         "equation carries no information; a plain Poisson GLM is the ",
         "appropriate model here.", call. = FALSE)
  }
  if (all(s == 0)) {
    stop("No unit is selected, so the outcome is never observed and the ",
         "outcome equation cannot be estimated.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check the Observed Counts
#'
#' @param y The response of the outcome equation.
#' @param selected Logical vector marking the selected units.
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_outcome_counts <- function(y, selected) {
  if (!is.numeric(y)) {
    stop("The outcome (left-hand side of 'outcome') must be a numeric count ",
         "variable, but is of class ", sQuote(class(y)[1L]), ".",
         call. = FALSE)
  }
  observed <- y[selected]
  if (anyNA(observed)) {
    stop("The outcome is missing for ", sum(is.na(observed)), " selected ",
         "unit(s). It must be observed wherever the selection indicator is 1; ",
         "NA is only allowed for non-selected units.", call. = FALSE)
  }
  is_invalid <- !is.finite(observed) | observed < 0 |
    observed != trunc(observed)
  if (any(is_invalid)) {
    stop("The outcome must consist of non-negative integer counts, but ",
         sum(is_invalid), " of the ", length(observed), " selected units ",
         "violate this.", call. = FALSE)
  }
  # If every observed count is 0 the ML estimate does not exist, the
  # intercept would just run off to -Inf. Better to stop here with a clear
  # message than to return an intercept of -27 with a huge standard error.
  if (all(observed == 0)) {
    stop("The outcome is 0 for every selected unit, so the outcome equation ",
         "cannot be estimated.", call. = FALSE)
  }
  report_ignored_outcomes(y, selected)
  invisible(TRUE)
}

#' Report Outcome Values that the Model Ignores
#'
#' Sometimes y is available for non-selected units too (simulated data often
#' has the complete outcome). That is not an error, the model just never uses
#' those values. We still print a message, because it could also mean that
#' the user mixed up the selection indicator.
#'
#' @inheritParams check_outcome_counts
#'
#' @return `invisible(TRUE)`.
#' @noRd
report_ignored_outcomes <- function(y, selected) {
  n_ignored <- sum(!selected & !is.na(y))
  if (n_ignored > 0L) {
    message("Note: the outcome is observed for ", n_ignored, " non-selected ",
            "unit(s). These values are ignored, because the model assumes ",
            "that the outcome is only observed when the selection indicator ",
            "is 1.")
  }
  invisible(TRUE)
}

#' Check a Design Matrix for Completeness and Full Rank
#'
#' @param design A numeric design matrix.
#' @param name The equation name, used in the error messages.
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_design_matrix <- function(design, name) {
  if (ncol(design) == 0L) {
    stop("The '", name, "' equation has no covariates. Please specify at ",
         "least an intercept.", call. = FALSE)
  }
  is_broken <- colSums(!is.finite(design)) > 0L
  if (any(is_broken)) {
    where <- switch(name, outcome = "for the selected units",
                    selection = "for all units")
    stop("The covariates of the '", name, "' equation must be complete ", where,
         " (no NA, NaN or Inf), but the following column(s) are not: ",
         toString(sQuote(colnames(design)[is_broken])), ".", call. = FALSE)
  }
  decomposition <- qr(design)
  if (decomposition$rank < ncol(design)) {
    # Nice trick: qr() with pivoting moves the linearly dependent columns to
    # the end of the pivot vector, so we can tell the user exactly which
    # columns are the problem instead of just saying "rank deficient".
    aliased <- decomposition$pivot[-seq_len(decomposition$rank)]
    stop("The design matrix of the '", name, "' equation is rank deficient, ",
         "so its coefficients are not identified. The following column(s) ",
         "are collinear with the remaining ones: ",
         toString(sQuote(colnames(design)[aliased])), ".", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check that the Sample is Large Enough
#'
#' @param model A model list as built by [build_model_data()].
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_sample_size <- function(model) {
  n_parameter <- model$n_beta + model$n_gamma + 2L
  if (model$n <= n_parameter) {
    stop("The model has ", n_parameter, " parameters but 'data' provides ",
         "only ", model$n, " observations.", call. = FALSE)
  }
  # beta and sigma only get information from the selected units, so we also
  # need enough of those, not just enough rows in total.
  if (model$n_selected <= model$n_beta + 1L) {
    stop("The outcome equation has ", model$n_beta, " coefficients plus ",
         "sigma, but the outcome is observed for only ", model$n_selected,
         " unit(s).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Warn if There is No Exclusion Restriction
#'
#' An exclusion restriction is a variable that affects the selection but is
#' not in the outcome equation (z1 in our simulated data). Strictly speaking
#' the model is identified without one, only through the non-linearity of
#' Phi, but in practice that is very weak: the profile of rho gets flat and
#' the standard errors explode. So we warn, but we do not stop.
#'
#' @param model A model list as built by [build_model_data()].
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_exclusion_restriction <- function(model) {
  # The intercept does not count: it shifts the selection probability for
  # everyone in the same way, so it cannot help to identify rho.
  exclusive <- setdiff(colnames(model$z), c("(Intercept)", colnames(model$x)))
  if (length(exclusive) == 0L) {
    warning("No exclusion restriction found: every covariate of the selection ",
            "equation also appears in the outcome equation. The model is then ",
            "only identified through the functional form, and the estimate ",
            "of 'rho' can be very unstable. Consider adding a variable that ",
            "affects the selection but not the outcome.", call. = FALSE)
  }
  invisible(TRUE)
}
