#' Reorder Data Blocks to Achieve Target Multivariate Correlation Structure
#'
#' This function reorders blocks of variables within a dataset so that the Spearman rank
#' correlations between the "representative" variables of each block match a specified target matrix.
#' It uses the Iman-Conover (1982) method. The function preserves the internal structure of each block
#' (intra-block correlations) and keeps the anchor block (and any unspecified columns) fixed in their original rows.
#'
#' @param data A data.frame containing the variables to be reordered.
#' @param blocks A named list defining the variable blocks.
#'   \itemize{
#'     \item **Names**: Must match the row/column names of \code{target_cor}.
#'     \item **Values**: Character vectors of column names in \code{data}.
#'     \item The **first variable** in each vector is treated as the "Representative" (Score).
#'           Its values determine the sorting order for that block and are used for correlation control.
#'     \item All variables in the vector move together as a unit to preserve intra-block relationships.
#'   }
#'   Example: \code{list(x = c("x", "x1"), m = c("m", "m1"))}.
#' @param target_cor A symmetric numeric matrix representing the desired Spearman rank correlation
#'   structure between the block representatives. Row and column names must match the names of the \code{blocks} list.
#' @param anchor A character string specifying which block name (from \code{blocks}) acts as the Anchor.
#'   The Anchor block is **never moved**.
#'
#' @return A data.frame with the same dimensions as \code{data}. Columns belonging to non-anchor blocks
#'   are reordered. Columns not specified in \code{blocks} (e.g., demographics) remain in their original positions.
#'
#' @importFrom MASS mvrnorm
#' @export
#'
#' @examples
#' # --- Example 1: Standard Usage ---
#'
#' # 1. Create dummy data with block structures
#' N <- 1000
#' df <- data.frame(
#'   # Demographics (Should not change relative to Anchor)
#'   id = 1:N,
#'   gender = sample(c("M", "F"), N, replace = TRUE),
#'   # Block X (Anchor): x is the score
#'   x = rnorm(N), x_item1 = rnorm(N), x_item2 = rnorm(N), x_item3 = rnorm(N),
#'   # Block M: m is the score
#'   m = runif(N), m_item1 = runif(N), m_item2 = runif(N), m_item3 = runif(N),
#'   # Block Y: y is the score
#'   y = rbeta(N, 2, 2), y_item1 = rbeta(N, 1, 1), y_item2 = rbeta(N, 1, 1), y_item3 = rbeta(N, 1, 1)
#' )
#' # 2. Define Block Structure
#' # Keys (x, m, y) match the target matrix.
#' # First value in vector is the representative.
#' block_def <- list(
#'   x = c("x", "x_item1", "x_item2", "x_item3"),
#'   m = c("m", "m_item1", "m_item2", "m_item3"),
#'   y = c("y", "y_item1", "y_item2", "y_item3")
#' )
#'
#' # 3. Define Target Matrix (x-m: 0.5, x-y: 0.3, m-y: 0.4)
#' target_mat <- matrix(c(
#'   1.0, 0.5, 0.3,
#'   0.5, 1.0, 0.4,
#'   0.3, 0.4, 1.0
#' ), nrow = 3, dimnames = list(c("x", "m", "y"), c("x", "m", "y")))
#'
#' # 4. Run function
#' df_new <- reorder_multivar(df, block_def, target_mat, anchor = "x")
#'
#' # 5. Validate
#' # Check correlations of representatives
#' print(cor(df_new[, c("x", "m", "y")], method = "spearman"))
#'
#' # Check if 'id' matches 'x' (Anchor didn't move)
#' print(all(df_new$x == df$x))     # TRUE
#' print(all(df_new$id == df$id))   # TRUE
#'
#' # Check intra-block integrity (m vs m_item1 correlation should be unchanged)
#' print(cor(df$m, df$m_item1))
#' print(cor(df_new$m, df_new$m_item1))
#'
#'
#' # --- Example 2: Non-Positive Semi-Definite Matrix (Error Handling) ---
#'
#' # Impossible structure: x strongly correlated to m and y (0.9),
#' # but m and y are negatively correlated (-0.9).
#' # This violates the triangle inequality.
#' bad_mat <- matrix(c(
#'   1.0, 0.9, 0.9,
#'   0.9, 1.0, -0.9,
#'   0.9, -0.9, 1.0
#' ), nrow = 3, dimnames = list(c("x", "m", "y"), c("x", "m", "y")))
#'
#'   # This will throw an error and provide a suggestion
#'   reorder_multivar(df, block_def, bad_mat, anchor = "x")
#'
reorder_multivar <- function(data, blocks, target_cor, anchor) {

  # --- 1. Input Validation ---

  if (!anchor %in% names(blocks)) {
    stop(sprintf("Anchor '%s' is not defined in the blocks list.", anchor))
  }

  if (!all(names(blocks) %in% colnames(target_cor))) {
    stop("Names in 'blocks' list must match row/col names in 'target_cor'.")
  }

  # Check if all columns exist in data
  all_cols <- unlist(blocks)
  missing_cols <- setdiff(all_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Columns not found in data:", paste(missing_cols, collapse = ", ")))
  }

  # Ensure target_cor is ordered correctly according to blocks list
  # This prevents mismatch if the user provided the matrix in a different order
  block_names <- names(blocks)
  target_cor <- target_cor[block_names, block_names]

  # --- 2. Positive Semi-Definite Check & Suggestion ---

  # Check eigenvalues to ensure the matrix is geometrically valid
  eigen_decomp <- eigen(target_cor, symmetric = TRUE)

  # Use a small tolerance for floating point comparisons
  if (any(eigen_decomp$values < -1e-8)) {

    # Attempt to find nearest valid correlation matrix by clipping negative eigenvalues
    # Method: Projection (Higham-like approach simplified)
    vals_clean <- pmax(eigen_decomp$values, 1e-7) # Floor at epsilon
    Q <- eigen_decomp$vectors
    mat_recon <- Q %*% diag(vals_clean) %*% t(Q)

    # Rescale to ensure diagonals are exactly 1 (convert covariance back to correlation)
    D_inv_sqrt <- diag(1 / sqrt(diag(mat_recon)))
    mat_nearest <- D_inv_sqrt %*% mat_recon %*% D_inv_sqrt

    # Preserve names
    dimnames(mat_nearest) <- dimnames(target_cor)

    # Format suggestion string
    suggest_str <- paste(utils::capture.output(print(round(mat_nearest, 3))), collapse = "\n")

    stop(paste0(
      "The target correlation matrix is not Positive Semi-Definite (it contains conflicting correlations).\n",
      "Mathematically, this structure cannot exist (e.g., violates triangle inequalities).\n\n",
      "SUGGESTION: The nearest mathematically valid matrix to your request is:\n",
      "------------------------------------------------------------------\n",
      suggest_str,
      "\n------------------------------------------------------------------\n",
      "Advice: You can copy the matrix above or reduce the magnitude of conflicting correlations."
    ))
  }

  # --- 3. Iman-Conover: Generate Reference Structure ---

  n_rows <- nrow(data)
  n_blocks <- length(blocks)

  # Generate multivariate normal samples with exact correlation structure
  # empirical = TRUE enforces exact correlation matching for the sample
  # We use MASS::mvrnorm directly as requested
  ref_data <- MASS::mvrnorm(n = n_rows, mu = rep(0, n_blocks), Sigma = target_cor, empirical = TRUE)
  colnames(ref_data) <- block_names

  # --- 4. Anchor Alignment ---

  # We cannot move the real anchor (data).
  # So, we must rearrange the *reference* matrix so that its anchor column
  # has the exact same rank order as the *real* anchor variable.

  # A. Get the representative variable name for the anchor block
  anchor_rep_var <- blocks[[anchor]][1]

  # B. Calculate Rank of the REAL anchor data
  # ties.method = "first" ensures a unique permutation vector 1..N
  real_anchor_rank <- rank(data[[anchor_rep_var]], ties.method = "first")

  # C. Sort the REFERENCE data by its anchor column
  # This creates a baseline where the reference anchor is strictly increasing
  ref_data_sorted <- ref_data[order(ref_data[, anchor]), , drop = FALSE]

  # D. Shuffle reference to match real anchor's rank structure
  # Now, ref_aligned[, anchor] is monotonic with data[[anchor_rep_var]]
  ref_aligned <- ref_data_sorted[real_anchor_rank, , drop = FALSE]

  # --- 5. Block Reordering ---

  result_df <- data # Start with original data (preserves id, gender, etc.)

  # Identify which blocks need moving (everything except the anchor)
  blocks_to_move <- setdiff(block_names, anchor)

  for (b_name in blocks_to_move) {

    # Get the columns belonging to this block
    cols_in_block <- blocks[[b_name]]
    rep_var <- cols_in_block[1] # The first one is the representative/score

    # Extract the chunk of real data for this block
    block_data <- data[, cols_in_block, drop = FALSE]

    # A. Sort the real block data by its representative variable
    #    This binds all sub-variables (e.g., m1, m2) to the score (m)
    #    ensuring they move as a cohesive unit.
    block_data_sorted <- block_data[order(block_data[[rep_var]]), , drop = FALSE]

    # B. Determine the target positions based on the aligned reference
    #    We look at the rank of the specific column in the reference matrix.
    #    This tells us: "For row i, we need the k-th smallest value from this block"
    target_rank <- rank(ref_aligned[, b_name], ties.method = "first")

    # C. Place the sorted real block into the positions dictated by the reference rank
    result_df[, cols_in_block] <- block_data_sorted[target_rank, , drop = FALSE]
  }

  return(result_df)
}
