#' Reorder Data Blocks in Place
#'
#' This function sorts the entire dataset based on an `anchor_score`, keeping
#' all row associations (IDs, unmentioned columns) intact. Then, for specific
#' `blocks` of columns, it shuffles them independently to achieve a target
#' correlation with the sorted anchor score.
#'
#' @param data A data.frame containing the original data.
#' @param anchor_score A string specifying the column name used to sort the entire dataset.
#'        All columns NOT in `blocks` will move together based on this sort.
#' @param blocks A named list of lists. Each sub-list defines a block to be reordered:
#'        \itemize{
#'          \item \code{cols}: Vector of column names for this block (items + total).
#'          \item \code{score}: The specific column used to determine internal order.
#'          \item \code{rho}: Target correlation with the \code{anchor_score} (-1 to 1).
#'        }
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return A data.frame with the same dimensions and column names as the input.
#'         Rows are sorted by `anchor_score`, and `blocks` columns are shuffled.
#' @export
reorder_block <- function(data, anchor_score, blocks, verbose = TRUE) {

  # --- Internal Helper: The Iman-Conover Logic ---
  # calculats indices to reorder 'block_score_vec' to correlate with 'target_vec'
  calc_new_indices <- function(target_vec, block_score_vec, rho) {
    n <- length(target_vec)
    # 1. Create structure for target
    target_rank <- rank(target_vec, ties.method = "random")
    target_signal <- scale(target_rank)

    # 2. Create noise structure
    noise <- rnorm(n)
    noise_signal <- scale(noise)

    # 3. Mix signal and noise
    # Formula: Target * rho + Noise * sqrt(1 - rho^2)
    mixed_signal <- target_signal * rho + noise_signal * sqrt(1 - rho^2)

    # 4. Determine mapping
    dest_order <- order(mixed_signal)   # Where the values should end up
    src_order  <- order(block_score_vec) # Where the values come from (sorted)

    indices <- integer(n)
    indices[dest_order] <- src_order
    return(indices)
  }

  # ============================================================
  # 1. Sort the ENTIRE Dataset (The "Anchor" Step)
  # ============================================================
  # This ensures ID, Gender, Y-block (if not processed), and X-block
  # all stay aligned row-wise according to the X score.
  if(verbose) message(sprintf("--- Step 1: Sorting entire dataset by anchor '%s' ---", anchor_score))

  # We work on a sorted copy of the data
  # This establishes the baseline: High Anchor Score -> Low Anchor Score
  ord <- order(data[[anchor_score]], decreasing = TRUE)
  data_sorted <- data[ord, ]

  # The target vector is now the sorted anchor score
  target_vector <- data_sorted[[anchor_score]]

  # ============================================================
  # 2. Process Blocks (The "Shuffle" Step)
  # ============================================================
  for(b_name in names(blocks)) {
    cfg <- blocks[[b_name]]

    if(verbose) {
      message(sprintf("--- Processing Block '%s' (Target rho = %.2f) ---", b_name, cfg$rho))
    }

    # Extract the block columns.
    # CRITICAL: We take these from the original data (or just treat them as a "bag of rows").
    # We need to find a permutation of these specific columns that matches the target_vector.
    raw_block <- data[, cfg$cols, drop = FALSE]

    # Calculate the new indices based on the Block's internal score
    new_idx <- calc_new_indices(
      target_vec = target_vector,          # The sorted X
      block_score_vec = raw_block[[cfg$score]], # The M score
      rho = cfg$rho                        # Desired correlation
    )

    # Reorder the block columns
    shuffled_block <- raw_block[new_idx, , drop = FALSE]

    # ============================================================
    # 3. Inject back into the sorted DataFrame
    # ============================================================
    # We overwrite the columns in data_sorted with the newly shuffled version
    data_sorted[, cfg$cols] <- shuffled_block
  }

  return(data_sorted)
}
