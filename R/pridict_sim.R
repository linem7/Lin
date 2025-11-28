#' Generate Predictors Based on an Existing Outcome Variable
#'
#' This function simulates categorical or ordinal predictor variables based on an existing
#' continuous or Likert-scale outcome variable. It uses a latent variable approach to ensure
#' the generated predictors maintain specific statistical relationships (correlation/significance)
#' with the outcome.
#'
#' @param y A numeric vector. The existing outcome variable (e.g., a Likert scale score).
#' @param predictor_configs A list of lists. Each sub-list represents a predictor to be generated
#'   and must contain:
#'   \itemize{
#'     \item \code{name}: String. The name of the new variable.
#'     \item \code{probs}: Numeric vector. The proportions for each category (e.g., \code{c(0.4, 0.6)}).
#'     \item \code{significant}: Logical. Whether the variable should be statistically significant.
#'     \item \code{effect_r}: (Optional) Numeric. The target correlation coefficient (-1 to 1).
#'       If provided, this takes precedence over the \code{significant} flag.
#'   }
#'
#' @return A data frame containing the original outcome variable and the newly generated predictors.
#' @importFrom stats rnorm residuals lm quantile sd runif as.formula
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' N <- 200
#' # Simulate an existing Likert outcome variable (1-7 scale)
#' existing_outcome <- sample(1:7, N, replace = TRUE,
#'                            prob = c(0.05, 0.1, 0.15, 0.3, 0.2, 0.15, 0.05))
#'
#' my_configs <- list(
#'   # Variable 1: Gender (2 categories, not significant but with slight specific noise)
#'   list(
#'     name = "Gender",
#'     probs = c(0.4, 0.6),      # Category 1: 40%, Category 2: 60%
#'     significant = FALSE,      # Not significant
#'     effect_r = 0.1            # Weak specific correlation
#'   ),
#'   # Variable 2: Education (5 categories, significant, default correlation)
#'   list(
#'     name = "Education",
#'     probs = rep(0.2, 5),      # Equal proportions
#'     significant = TRUE,       # Significant
#'     effect_r = 0.4            # Medium correlation
#'   ),
#'   # Variable 3: Income Group (3 categories, custom props, strong correlation)
#'   list(
#'     name = "Income_Group",
#'     probs = c(0.2, 0.5, 0.3), # 20%, 50%, 30%
#'     significant = TRUE,
#'     effect_r = 0.65           # Strong correlation
#'   )
#' )
#'
#' final_data <- generate_predictors(existing_outcome, my_configs)
#' validate_simulation(final_data, outcome_col = "outcome")
#' }
generate_predictors <- function(y, predictor_configs) {

  # Check input
  if (!is.numeric(y)) stop("Outcome variable y must be numeric (Likert)")
  n <- length(y)

  # Initialize result data frame
  df_sim <- data.frame(outcome = y)

  # Scale y for calculation (does not affect original values)
  y_scaled <- scale(y)

  # Loop through configs to generate each variable
  for (conf in predictor_configs) {

    # --- Step 1: Determine effect size (Correlation) ---

    # Initialize target_r
    target_r <- 0

    # Logic: Prioritize effect_r if specified
    if (!is.null(conf$effect_r)) {
      # Case A: User specified strength -> Strictly follow, ignore significant flag
      target_r <- conf$effect_r

    } else {
      # Case B: User did not specify strength -> Auto-generate based on significance
      if (conf$significant) {
        # Auto-generate significant correlation (default 0.3 ~ 0.6)
        target_r <- stats::runif(1, 0.3, 0.6) * sample(c(1, -1), 1)
      } else {
        # Auto-generate non-significant noise (default -0.05 ~ 0.05)
        target_r <- stats::runif(1, -0.05, 0.05)
      }
    }

    # --- Step 2: Generate Latent Variable ---
    # Formula: X_latent = r * Y + sqrt(1 - r^2) * Noise
    # This ensures X_latent has approx correlation 'target_r' with Y
    noise <- stats::rnorm(n)

    # Orthogonalization: Ensure noise and y are truly uncorrelated for mathematical precision
    resid_noise <- stats::residuals(stats::lm(noise ~ y_scaled))
    resid_noise <- scale(resid_noise) # Re-scale

    x_latent <- target_r * y_scaled + sqrt(1 - target_r^2) * resid_noise

    # --- Step 3: Discretize based on proportions ---
    # Normalize proportions
    probs <- conf$probs / sum(conf$probs)

    # Calculate cumulative probabilities for cut points
    cut_points <- stats::quantile(x_latent, probs = cumsum(probs), names = FALSE)

    # Add -Inf and cut points (remove the last 100% point to prevent boundary issues)
    breaks <- c(-Inf, cut_points[-length(cut_points)], Inf)

    # Cut and convert to integer (1, 2, 3...)
    x_discrete <- cut(x_latent, breaks = breaks, labels = FALSE)

    # Store in data frame
    df_sim[[conf$name]] <- x_discrete
  }

  return(df_sim)
}

#' Validate Simulated Data
#'
#' Runs a linear regression model and checks frequency distributions to validate
#' the generated predictor variables against the outcome.
#'
#' @param df A data frame containing the outcome and predictors.
#' @param outcome_col A string specifying the name of the outcome column. Defaults to "outcome".
#'
#' @return A list containing the linear model object and the data frame.
#' @importFrom stats lm as.formula summary.lm
#' @export
validate_simulation <- function(df, outcome_col = "outcome") {
  cat("================ Simulation Validation Report ================\n\n")

  # 1. Run regression model
  predictors <- setdiff(names(df), outcome_col)
  f <- stats::as.formula(paste(outcome_col, "~", paste(predictors, collapse = " + ")))
  model <- stats::lm(f, data = df)

  model_sum <- summary(model)

  cat("--- Regression Analysis Results (Linear Regression) ---\n")
  print(model_sum$coefficients)
  cat("\nR-squared:", round(model_sum$r.squared, 3), "\n\n")

  # 2. Validate category distributions
  cat("--- Category Proportion Validation ---\n")
  for (p in predictors) {
    cat(sprintf("Variable [%s] Distribution:\n", p))
    tbl <- table(df[[p]])
    props <- prop.table(tbl)
    print(round(props, 3))
    cat("\n")
  }

  return(list(model = model, data = df))
}
