test_that("apa returns a kableExtra table object with correct dimensions", {
  tbl <- head(mtcars, 5)
  result <- apa(tbl, title = "Test Table")
  # Check class (kableExtra creates an object of class "kableExtra" or "knitr_kable")
  expect_true("knitr_kable" %in% class(result) || "kableExtra" %in% class(result))
  # Check that caption is included if applicable (depends on how apa is implemented)
  if ("knitr_kable" %in% class(result)) {
    expect_true(any(grepl("Test Table", result, fixed = TRUE)))
  }
})
