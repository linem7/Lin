#' Plot mixture‐model item‐response profiles
#'
#' Extracts conditional means and/or response probabilities from one or more
#' MplusAutomation model objects and produces a ggplot2 profile plot across
#' latent profiles.
#'
#' @param models A single MplusAutomation model object or a list of such objects.
#' @param model_labels Optional character vector of labels for each model (e.g.
#'   "3 Classes"). If NULL, labels will be auto‐generated from the number of
#'   classes in each model.
#' @param select Integer indices or names of models to keep. NULL (default)
#'   retains all models.
#' @param param_order Optional character vector specifying the desired order
#'   of the parameters on the x‐axis. If NULL and `param_rename` is provided,
#'   the new names will be used in the order they appear in `param_rename`.
#' @param param_rename Optional named character vector for recoding parameter
#'   names, e.g. `c(oldName = "New Name")`. The order of values here will be
#'   used as the default `param_order` if not otherwise specified.
#' @param prob_category Integer specifying which category to extract from
#'   `parameters$probability.scale` (defaults to 2).
#' @param param_exclude Optional character vector of parameter names to drop
#'   (e.g. `c("DAYS_GAP")`).
#' @param return_data Logical; if TRUE returns the processed data frame instead
#'   of plotting (default FALSE).
#'
#' @return A `ggplot` object showing item‐response profiles across latent
#'   profiles, or (if `return_data = TRUE`) a data frame with columns:
#'   `Model`, `Profiles`, `label`, `param`, and `est`.
#'
#' @import dplyr
#' @import purrr
#' @import ggplot2
#' @examples
#' # 1. Load the pre‐read model objects you saved in data/lpa_res.RData
#' data(lpa_res)
#'
#' # 1. Multi‐model plot (auto‐generated "3 Classes", "4 Classes", etc.)
#' plot_mixture_profile(models = lpa_res)
#'
#' # 2. Reorder parameters on the x‐axis to X9–X1
#' plot_mixture_profile(
#'   models      = lpa_res,
#'   param_order = paste0("X", 9:1)
#' )
#'
#' # 3. Rename X1–X9 to more descriptive labels
#' rename_map <- setNames(
#'   paste("Indicator", 1:9),
#'   paste0("X", 1:9)
#' )
#' plot_mixture_profile(
#'   models       = lpa_res,
#'   param_rename = rename_map
#' )
#'
#' # 4. Single‐model plot: legend will show "Class 1 (xx.xx%)", etc.
#' plot_mixture_profile(models = lpa_res$lpa_3c.out)
#'
#' # 5. Return data rather than a plot object
#' plot_mixture_profile(models = lpa_res$lpa_3c.out, return_data = T)
#'
#' @export
plot_mixture_profile <- function(
    models,
    model_labels  = NULL,
    select        = NULL,
    param_order   = NULL,
    param_rename  = NULL,
    prob_category = 2,
    param_exclude = NULL,
    return_data   = FALSE
) {
  # If single model object passed directly, wrap in list
  if (!is.null(models[["parameters"]]) && !is.null(models[["class_counts"]])) {
    models <- list(models)
  }

  n_models <- length(models)
  # Auto‐generate labels if needed
  if (is.null(model_labels)) {
    model_labels <- purrr::map_chr(models, ~ {
      nclass <- length(unique(.x[["class_counts"]][["mostLikely"]]$class))
      paste0(nclass, " Classes")
    })
  }
  stopifnot(length(model_labels) == n_models)

  # Helper: extract one model's data
  extract_one <- function(model, model_label) {
    # A. Means: drop latent‐variable means if none remain
    tmp_means <- model[["parameters"]][["unstandardized"]] %>%
      dplyr::filter(
        paramHeader == "Means",
        LatentClass != "Categorical.Latent.Variables"
      ) %>%
      dplyr::select(param, est, LatentClass)
    means_df <- if (nrow(tmp_means) > 0) tmp_means else NULL

    # B. Probabilities: only if table exists
    prob_df <- if ("probability.scale" %in% names(model[["parameters"]])) {
      model[["parameters"]][["probability.scale"]] %>%
        dplyr::filter(category == prob_category) %>%
        dplyr::select(param, est, LatentClass)
    } else {
      NULL
    }

    # C. Combine available data
    df <- dplyr::bind_rows(means_df, prob_df)

    # D. Build profile labels
    labels <- model[["class_counts"]][["mostLikely"]] %>%
      dplyr::transmute(
        LatentClass = as.character(class),
        label       = paste0(
          "Class ", class, " (",
          sprintf("%0.2f%%", proportion * 100),
          ")"
        )
      )

    # E. Join, rename and select final columns
    df %>%
      dplyr::left_join(labels, by = "LatentClass") %>%
      dplyr::mutate(
        Model    = model_label,
        Profiles = LatentClass  # rename for plotting legend
      ) %>%
      dplyr::select(Model, Profiles, label, param, est)
  }

  # Combine all models' data
  all_data <- purrr::map2(models, model_labels, extract_one) %>%
    dplyr::bind_rows()

  # Exclude specified parameters
  if (!is.null(param_exclude)) {
    all_data <- dplyr::filter(all_data, ! param %in% param_exclude)
  }
  # Subset models if requested
  if (!is.null(select)) {
    keep <- if (is.numeric(select)) model_labels[select] else select
    all_data <- dplyr::filter(all_data, Model %in% keep)
  }

  # Recode and reorder parameters
  if (!is.null(param_rename)) {
    all_data <- dplyr::mutate(all_data,
                              param = dplyr::recode(param, !!!param_rename))
  }

  # If param_order is NULL but user provided param_rename, derive order
  if (is.null(param_order) && !is.null(param_rename)) {
    param_order <- unname(param_rename)
  }
  if (!is.null(param_order)) {
    all_data <- dplyr::mutate(all_data,
                              param = factor(param, levels = param_order))
  }

  # Return data frame if requested
  if (return_data) {
    return(all_data)
  }

  # Build profile plot (multi‐ and single‐model uses same structure)
  multi <- length(unique(all_data$Model)) > 1
  if (multi) {
    p <- ggplot2::ggplot(
      all_data,
      ggplot2::aes(x = param, y = est, group = Profiles),
      color = "black"
    ) +
      ggplot2::geom_line(ggplot2::aes(lty = Profiles), color = "black") +
      ggplot2::geom_point(ggplot2::aes(shape = Profiles), color = "black") +
      ggplot2::facet_wrap(~ Model) +
      ggplot2::labs(
        title = "Figure x. Plot of conditional item response means.",
        x     = "Indicators",
        y     = "Estimates",
        linetype = "Profiles",
        shape    = "Profiles"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)
      )
  } else {
    p <- ggplot2::ggplot(
      all_data,
      ggplot2::aes(
        x = param, y = est,
        group = label, color = label, linetype = label, shape = label
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::labs(
        title    = "Figure x. Plot of conditional item response means.",
        x        = "Indicators",
        y        = "Estimates",
        color    = "Profiles",
        linetype = "Profiles",
        shape    = "Profiles"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5)
      )
  }
  p
}
