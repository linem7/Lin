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
#'
#' @examples
#'
#' create_cfa_model("pasu_f", "pasu{i}_f", 5)
#'
#' # For multiple factors, combine them with paste after place into a vector:
#' paste(c(create_cfa_model("pasu_f", "pasu{i}_f", 5),
#' create_cfa_model("cmu_f", "cmu{i}_f", 3),
#' create_cfa_model("att_f", "att{i}_f", 4)), collapse = "\n")
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
