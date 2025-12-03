library(tidyverse)
library(nnet) # Required for multinomial logistic regression

# ==============================================================================
# Core Function 1: Generate Predictors for Logistic/Multinomial Regression
# ==============================================================================

#' Generate Predictors for Logistic or Multinomial Regression
#'
#' This function simulates predictor variables (continuous, factor, or binary) based on an
#' existing categorical outcome variable (binary or multinomial). It uses a latent variable
#' approach where the mean of the latent variable differs across outcome groups to create
#' statistical associations.
#'
#' @param y A vector (numeric or factor). The existing categorical outcome variable.
#' @param predictor_configs A list of lists. Each sub-list defines a predictor to be generated
#'   and must contain:
#'   \itemize{
#'     \item \code{name}: String. The name of the new variable.
#'     \item \code{probs}: Numeric vector. The global proportions for each category of the predictor.
#'       These proportions are maintained across the entire dataset.
#'     \item \code{significant}: Logical. Whether the variable should be statistically significant
#'       (i.e., different latent means across outcome groups).
#'     \item \code{type}: String. The type of the variable in the model:
#'       \code{"continuous"} (numeric), \code{"factor"} (categorical), or \code{"binary"}.
#'     \item \code{group_means}: (Optional) Numeric vector. Specifies the mean of the latent variable
#'       for each category of the outcome \code{y}. The length must match the number of outcome levels.
#'       If provided, this overrides the \code{significant} flag.
#'   }
#'
#' @return A data frame containing the original outcome variable (as a factor) and the
#'   newly generated predictors.
#' @importFrom stats rnorm quantile runif as.formula glm binomial pnorm model.matrix
#' @importFrom nnet multinom
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(2023)
#' N <- 500
#' # Outcome variable Y with 3 categories (e.g., Brand A, B, C)
#' y_nominal <- sample(1:3, N, replace = TRUE, prob = c(0.3, 0.4, 0.3))
#'
#' my_logistic_configs <- list(
#'
#'   # Variable 1: Satisfaction Score (5 categories)
#'   # Treated as "continuous" in the model (linear trend 1-5)
#'   list(
#'     name = "Satisfaction_Score",
#'     probs = c(0.1, 0.2, 0.4, 0.2, 0.1),
#'     significant = TRUE,
#'     type = "continuous",
#'     # Linear trend: Group 1 low, Group 3 high
#'     group_means = c(-1, 0, 1)
#'   ),
#'
#'   # Variable 2: Region (3 categories)
#'   # Treated as "factor" (Dummy variables)
#'   list(
#'     name = "Region",
#'     probs = c(0.3, 0.3, 0.4),
#'     significant = TRUE,
#'     type = "factor",
#'     # Distinct distribution for each group
#'     group_means = c(0, 0.8, -0.5)
#'   ),
#'
#'   # Variable 3: Gender (2 categories)
#'   # Not significant, but with very minor manual differences
#'   list(
#'     name = "Gender",
#'     probs = c(0.5, 0.5),
#'     significant = FALSE,
#'     type = "factor",
#'     # Tiny differences to avoid mathematical equality (simulation realism)
#'     group_means = c(0, 0.1, -0.05)
#'   )
#' )
#'
#' # 1. Generate data
#' sim_data <- generate_logistic_predictors(y_nominal, my_logistic_configs)
#'
#' # 2. Validate data
#' validate_logistic_simulation(sim_data)
#' }
generate_logistic_predictors <- function(y, predictor_configs) {

  # Force outcome to factor (Grouping Variable)
  y_fac <- as.factor(y)
  levels_y <- levels(y_fac)
  n_groups <- length(levels_y)
  n <- length(y)

  # Initialize result data frame
  df_sim <- data.frame(outcome = y_fac)

  for (conf in predictor_configs) {

    # --- Step 1: Set Group Differences (Effect Size) ---
    # Logic: If X predicts Y, the mean of the latent variable for X differs across Y groups.

    group_means <- rep(0, n_groups)

    # Logic check: Prioritize user-provided group_means
    if (!is.null(conf$group_means)) {
      # Case A: User specified means -> Strictly follow, ignore 'significant' flag
      if(length(conf$group_means) != n_groups) {
        stop(paste0("Variable ", conf$name, ": length of group_means must equal number of Y categories"))
      }
      group_means <- conf$group_means

    } else {
      # Case B: User did not specify -> Auto-generate based on 'significant' flag
      if (conf$significant) {
        # Auto-generate significant differences: Random means between -0.8 and 0.8
        group_means <- stats::runif(n_groups, -0.8, 0.8)
        # Optionally force one group to be very distinct to ensure significance
        group_means[sample(1:n_groups, 1)] <- 1.5
      } else {
        # significant = FALSE
        # Generate tiny random noise (e.g., -0.05 to 0.05) to avoid exact mathematical equality
        group_means <- stats::runif(n_groups, -0.05, 0.05)
      }
    }

    # --- Step 2: Generate Latent Variable ---
    # Assign the group mean to each observation based on its Y value
    obs_means <- group_means[as.numeric(y_fac)]

    # Generate normal latent variable: Mean depends on Y, SD is fixed at 1
    x_latent <- stats::rnorm(n, mean = obs_means, sd = 1)

    # --- Step 3: Discretize based on Global Proportions ---
    # We must maintain the proportions specified in conf$probs across the WHOLE dataset
    probs <- conf$probs / sum(conf$probs)

    # Calculate quantiles based on the entire x_latent distribution
    cut_points <- stats::quantile(x_latent, probs = cumsum(probs), names = FALSE)

    # Handle boundaries and cut
    breaks <- c(-Inf, cut_points[-length(cut_points)], Inf)
    x_discrete <- cut(x_latent, breaks = breaks, labels = FALSE)

    # --- Step 4: Convert variable type as requested ---
    if (conf$type == "continuous") {
      # Treated as continuous numeric (1, 2, 3...)
      df_sim[[conf$name]] <- as.numeric(x_discrete)

    } else if (conf$type == "factor") {
      # Treated as categorical factor
      df_sim[[conf$name]] <- as.factor(x_discrete)

    } else if (conf$type == "binary") {
      # Binary variable: usually 0/1 or 1/2. Here we return as Factor.
      df_sim[[conf$name]] <- as.factor(x_discrete)
    } else {
      # Default to numeric if unspecified
      df_sim[[conf$name]] <- as.numeric(x_discrete)
    }
  }

  return(df_sim)
}


# ==============================================================================
# Core Function 2: Validation (Supports Binary and Multinomial Logistic Regression)
# ==============================================================================

#' Validate Logistic Simulation
#'
#' fits a logistic regression (binary) or multinomial logistic regression model
#' to the simulated data and prints a summary of coefficients, standard errors,
#' and p-values to verify the generated relationships.
#'
#' @param df A data frame containing the outcome and predictors.
#' @param outcome_col String. The name of the outcome column. Defaults to "outcome".
#'
#' @return NULL. The function prints the validation report to the console.
#' @export
validate_logistic_simulation <- function(df, outcome_col = "outcome") {
  cat("================ Logistic Regression Validation Report ================\n\n")

  y <- df[[outcome_col]]
  predictors <- setdiff(names(df), outcome_col)
  f <- stats::as.formula(paste(outcome_col, "~", paste(predictors, collapse = " + ")))

  n_levels <- length(unique(y))

  # --- 1. Run Regression Model ---
  if (n_levels == 2) {
    cat("Detected binary outcome variable, using Binary Logistic Regression (glm-binomial)...\n")
    model <- stats::glm(f, data = df, family = binomial)
    print(summary(model))

  } else {
    cat("Detected multinomial outcome variable (", n_levels, " levels), using Multinomial Logistic Regression (multinom)...\n")
    # trace=FALSE hides iteration progress
    model <- nnet::multinom(f, data = df, trace = FALSE)

    sum_mod <- summary(model)
    coefs <- sum_mod$coefficients
    errors <- sum_mod$standard.errors

    # multinom summary returns matrices (rows=contrasts, cols=predictors)
    # We need to generate a summary table similar to summary(glm) for each contrast group
    contrast_groups <- rownames(coefs)

    for (i in 1:nrow(coefs)) {
      group_name <- contrast_groups[i]
      cat(sprintf("\n=== Comparison: %s vs Reference (1) ===\n", group_name))

      est <- coefs[i, ]
      se <- errors[i, ]
      z_val <- est / se
      p_val <- (1 - stats::pnorm(abs(z_val))) * 2

      # Build result matrix
      res_matrix <- cbind(
        Estimate = est,
        `Std. Error` = se,
        `z value` = z_val,
        `Pr(>|z|)` = p_val
      )

      # Use printCoefmat to automatically format and add significance stars
      stats::printCoefmat(res_matrix, P.values = TRUE, has.Pvalue = TRUE, signif.stars = TRUE)
    }
  }

  # --- 2. Validate Predictor Proportions ---
  cat("\n================ Variable Frequency Validation ================\n")
  cat(sprintf("Total Observations (N): %d\n", nrow(df))) # Added Total N for reference

  for (p in predictors) {
    cat(sprintf("\nDistribution of variable [%s] (Type: %s):\n", p, class(df[[p]])))

    # Generate the frequency table (Counts)
    tbl <- table(df[[p]])

    # Print the table directly to show counts instead of proportions
    print(tbl)
  }
}
