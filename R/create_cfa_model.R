#' Create lavaan model
#'
#' A more easy way to create lavaan style cfa model syntax
#'
#' @param lav_name names of latent variable
#' @param pattern pattern of the model structure
#' @param var_index which item should be counted into the model
#'
#' @returns a string contains the model structure
#' @export
#' @importFrom glue glue glue_collapse
#' @examples
#' # 1) Single-factor example:
#' model1 <- create_cfa_model("x", "x{i}", 1:3)
#' # which would return model1 <- 'x =~ x1 + x2 + x3`
#'
#' # Evaluate fit and reliability:
#' lavaan::cfa(model1, data = data)
#' reliability_table(model1, data = mydata, digits = 3)
#'
#' # 2) Multi-factor model, then reliability:
#' model2 <- paste(
#'   create_cfa_model("x", "x{i}", 1:5),
#'   create_cfa_model("m", "m{i}", 1:3),
#'   create_cfa_model("z", "z{i}", 1:4),
#'   collapse = "\n"
#' )
#' # Now get reliability for the combined model:
#' reliability_table(model2, data = mydata, digits = 3)
#'
#'
create_cfa_model <- function(lav_name, pattern, var_index) {
  # var_index can be a vector like c(1,7,8)
  # Generate items using pattern and specified indices
  items <- sapply(var_index, function(i) {
    glue(pattern, i = i)
  })

  # Create model syntax
  model <- glue("{lav_name} =~ {glue_collapse(items, ' + ')}")

  return(model)
}
