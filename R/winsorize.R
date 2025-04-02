#' Helper function to recode outliers
#'
#' Winsorizing value which outside the 3 SD boundary.
#' @param x column of a data.frame
#'
#' @returns a new columns after mutate function
#' @export
#'
#' @examples
#' data_clean <- data %>%
#' mutate(
#' # Recode rt and acc value
#'       rt_recode = winsorize(rt_m),
#'       acc_recode = winsorize(acc_m))
#'
#'
#'
winsorize <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  upper_limit <- mean_x + 3 * sd_x
  lower_limit <- mean_x - 3 * sd_x

  case_when(
    x > upper_limit ~ upper_limit,
    x < lower_limit ~ lower_limit,
    TRUE ~ x
  )
}
