#' Reorder Data Blocks to Simulate Target Correlations
#'
#' This function reorders blocks of variables (e.g., items + total score) to match
#' a specified target correlation with an "Anchor" variable. It moves the entire
#' block (items and scores) together to preserve internal consistency (Cronbach's alpha).
#'
#' @param data A data.frame containing the original data.
#' @param anchor_cols A character vector of column names that form the "Anchor" block.
#'        These variables are sorted first and define the reference order.
#'        (e.g., demographics + independent variable items + total score).
#' @param anchor_score A string specifying the column name within \code{anchor_cols}
#'        used to sort the anchor block.
#' @param blocks A named list of lists. Each sub-list defines a block to be reordered:
#'        \itemize{
#'          \item \code{cols}: Vector of column names for this block (items + total).
#'          \item \code{score}: The specific column used to determine internal order.
#'          \item \code{rho}: Target correlation with the \code{anchor_score} (-1 to 1).
#'        }
#' @param verbose Logical. If \code{TRUE}, prints validation checks.
#'
#' @return A data.frame with shuffled rows where correlations match the targets.
#' @export
#'
#' @examples
#' # --- 1. Create Mock Data ---
#' set.seed(123)
#' N <- 1000
#' d <- data.frame(
#'   # Continuous Items (x1-x3, m1-m3, y1-y3)
#'   x1 = rnorm(N), x2 = rnorm(N), x3 = rnorm(N),
#'   m1 = rnorm(N), m2 = rnorm(N), m3 = rnorm(N),
#'   y1 = rnorm(N), y2 = rnorm(N), y3 = rnorm(N),
#'
#'   # Demographics / Categorical
#'   gender = sample(c(0, 1), N, replace = TRUE),
#'   grade  = sample(c("k1", "k2", "k3"), N, replace = TRUE)
#' )
#'
#' # Calculate Total Scores
#' d$x <- rowSums(d[, c("x1","x2","x3")])
#' d$m <- rowSums(d[, c("m1","m2","m3")])
#' d$y <- rowSums(d[, c("y1","y2","y3")])
#'
#' # --- 2. Define Configuration ---
#' cols_anchor <- c("x1", "x2", "x3", "x", "gender", "grade")
#' block_config <- list(
#'   mediator_block = list(
#'     cols  = c("m1", "m2", "m3", "m"),
#'     score = "m",
#'     rho   = 0.60
#'   ),
#'   outcome_block = list(
#'     cols  = c("y1", "y2", "y3", "y"),
#'     score = "y",
#'     rho   = 0.30
#'   )
#' )
#'
#' # --- 3. Execute ---
#' data_final <- reorder_block(
#'   data = d,
#'   anchor_cols = cols_anchor,
#'   anchor_score = "x",
#'   blocks = block_config,
#'   verbose = TRUE
#' )
#'
#' # --- 4. Verify ---
#' # Check if correlations match the targets (approx 0.60 and 0.30)
#' cor(data_final[, c("x", "m", "y")])
reorder_block <- function(data, anchor_cols, anchor_score, blocks, verbose = TRUE) {

  # --- Internal Helper: The Iman-Conover Logic ---
  calc_new_indices <- function(target_vec, block_score_vec, rho) {
    n <- length(target_vec)
    target_rank <- rank(target_vec, ties.method = "random")
    target_signal <- scale(target_rank)
    noise <- rnorm(n)
    noise_signal <- scale(noise)
    mixed_signal <- target_signal * rho + noise_signal * sqrt(1 - rho^2)
    dest_order <- order(mixed_signal)
    src_order  <- order(block_score_vec)
    indices <- integer(n)
    indices[dest_order] <- src_order
    return(indices)
  }

  # ============================================================
  # 1. Anchor Step
  # ============================================================
  if(verbose) message("--- Step 1: Processing Anchor Block ---")

  anchor_order <- order(data[[anchor_score]], decreasing = TRUE)
  df_anchor <- data[anchor_order, anchor_cols, drop = FALSE]
  target_vector <- df_anchor[[anchor_score]]

  final_pieces <- list(anchor = df_anchor)

  # ============================================================
  # 2. Block Reordering Step
  # ============================================================
  for(b_name in names(blocks)) {
    cfg <- blocks[[b_name]]

    if(verbose) {
      message(sprintf("--- Processing %s (Target rho = %.2f) ---", b_name, cfg$rho))
    }

    raw_block <- data[, cfg$cols, drop = FALSE]

    new_idx <- calc_new_indices(
      target_vec = target_vector,
      block_score_vec = raw_block[[cfg$score]],
      rho = cfg$rho
    )

    reordered_block <- raw_block[new_idx, , drop = FALSE]
    final_pieces[[b_name]] <- reordered_block
  }

  # ============================================================
  # 3. Merge (Fixed to prevent renaming columns)
  # ============================================================
  # Unname the list so cbind doesn't prefix columns with block names
  result_df <- do.call(cbind, unname(final_pieces))

  return(result_df)
}
