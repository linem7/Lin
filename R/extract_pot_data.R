#' Exact plot data
#'
#' Exact plot data which combines both means and probabilities
#'
#' @param model object returned by \code{readModels} for mixture model
#' @param model_label index of the latent class model, e.g., 1-class model, 2-class model
#'
#' @returns a data.frame
#' @export
#'
#' @examples
#'
#'
#' # Applying the function to all model objects within the list
#' plot_data_pre <- map2(res_pre, 1:5, extract_plot_data) %>%
#'   bind_rows() %>%
#'   filter(!param == "DAYS_GAP") %>%
#'   mutate(param = recode(param,
#'                         "TEC_RECOG" = "Recognition",
#'                         "TEC_EXT" = "External causes",
#'                         "TEC_DES" = "Desire",
#'                         "TEC_BEL" = "Belief",
#'                         "TEC_REM" = "Reminder",
#'                         "TEC_REG" = "Regulation",
#'                         "TEC_HID" = "Hiding"),
#'          param = factor(param, levels = c("Recognition", "External causes", "Desire", "Belief", "Reminder", "Regulation", "Hiding"))) # Change the order of variables for plotting

extract_plot_data <- function(model, model_label) {
  # Extracting 'Means' data
  # The results stored in different place when using createMixtures
  # use model[["parameters"]][["unstandardized"]] when the result imported by readModels
  means_data <- model[["parameters"]][["unstandardized"]] %>%
    filter(paramHeader == "Means", !str_detect(LatentClass, "^Categorical")) %>%
    select(param, est, LatentClass)

  # Extracting 'probability.scale' data
  prob_scale_data <- model[["parameters"]][["probability.scale"]] %>%
    filter(category == 2) %>%
    select(param, est, LatentClass)

  # Extract class counts and proportions
  class_proportions <- model[["class_counts"]][["mostLikely"]] %>%
    mutate(proportion = str_c(formatC(proportion * 100, format = "f", digits = 2), "%")) %>%
    mutate(label = str_c("Class ", class, " (", proportion, ")"),
           LatentClass = as.character(class)) %>%
    select(LatentClass, label)

  # Join the proportions data with the means and prob data
  means_data <- means_data %>%
    left_join(class_proportions, by = "LatentClass")
  prob_scale_data <- prob_scale_data %>%
    left_join(class_proportions, by = "LatentClass")

  # Combine the data frames
  combined_data <- bind_rows(means_data, prob_scale_data) %>%
    mutate(Model = paste(model_label, "Classes")) %>%
    select(Model, LatentClass, label, param, est)

  return(combined_data)
}
