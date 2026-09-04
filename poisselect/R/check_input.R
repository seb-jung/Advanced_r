# Input checks, in two stages ("fail fast, fail appropriately"):
# check_arguments() validates what the user typed before any data is touched,
# check_model_data() validates what the formulas imply once the design matrices
# exist, which is where the NA pattern, collinearity and the exclusion
# restriction become visible.

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
  # A single node makes the rate independent of sigma, so K = 1 is rejected;
  # beyond a few dozen nodes the accuracy gain is nil while the cost grows.
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
  # '.' would put the response of the other equation among the covariates.
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
#' Runs after [build_model_data()], because the required length is only known
#' once the design matrices exist.
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
  # The outcome covariates enter the likelihood for the selected units only,
  # the selection covariates for every unit.
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
  # With only zeros the intercept of the outcome equation runs off to -Inf.
  if (all(observed == 0)) {
    stop("The outcome is 0 for every selected unit, so the outcome equation ",
         "cannot be estimated.", call. = FALSE)
  }
  report_ignored_outcomes(y, selected)
  invisible(TRUE)
}

#' Report Outcome Values that the Model Ignores
#'
#' An observed outcome of a non-selected unit is not an error, simulated data
#' regularly carries the complete outcome, but the value does not enter the
#' likelihood, so the user is told.
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
    # qr() pivots the linearly dependent columns to the end, which is how they
    # can be named in the message.
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
  # beta and sigma are estimated from the selected units alone.
  if (model$n_selected <= model$n_beta + 1L) {
    stop("The outcome equation has ", model$n_beta, " coefficients plus ",
         "sigma, but the outcome is observed for only ", model$n_selected,
         " unit(s).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Warn if There is No Exclusion Restriction
#'
#' Formally the model is identified through the non-linearity of the normal
#' distribution function alone, but that identification is weak. A variable
#' that shifts the selection probability without entering the outcome equation
#' is what makes the estimate of `rho` trustworthy in practice.
#'
#' @param model A model list as built by [build_model_data()].
#'
#' @return `invisible(TRUE)`.
#' @noRd
check_exclusion_restriction <- function(model) {
  # The intercept does not count as an exclusion restriction.
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
