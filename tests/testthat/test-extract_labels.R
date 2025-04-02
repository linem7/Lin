test_that("extract_labels returns correct labels for each column", {
  # Create a sample data frame with labels
  df <- data.frame(col1 = 1:3, col2 = 4:6)
  attr(df$col1, "label") <- "Column 1 Label"
  attr(df$col2, "label") <- "Column 2 Label"

  # Run the function
  result <- extract_labels(df)

  # Check that result is a data frame or tibble with the expected content
  expect_true(is.data.frame(result) || tibble::is_tibble(result))
  expect_equal(nrow(result), 2)                   # two rows for two columns
  expect_equal(as.character(result$variable), c("col1", "col2"))
  expect_equal(as.character(result$label), c("Column 1 Label", "Column 2 Label"))
})
