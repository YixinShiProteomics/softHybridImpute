#' Soft Hybrid Imputation for Proteomics Missing Values
#'
#' @description
#' Blends Random Forest (RF) and MinProb imputations with a soft weight driven by
#' missing rate (MNAR proxy) and mean intensity (abundance). The missing-rate elbow
#' is estimated by LOESS + farthest-point (Kneedle-like) unless 75\% of proteins
#' have zero missing (fallback to 1/n rule).
#'
#' @param df_rf A matrix/data.frame of RF-imputed values (rows = features, cols = samples).
#' @param df_minProb A matrix/data.frame of MinProb-imputed values (same dim/rownames as \code{df_rf}).
#' @param df_raw The raw matrix/data.frame with NAs (used to compute missing rate & mean intensity).
#' @param r0 Optional numeric, elbow x (missing-rate) override. If NULL, auto-estimated.
#' @param x0 Optional numeric, elbow y (mean-intensity) override. If NULL, auto-estimated.
#' @param a Steepness for missing-rate sigmoid (default 10).
#' @param b Steepness for intensity sigmoid (default 5).
#' @param lambda Trade-off between missing-rate and intensity (0–1, default 0.5).
#' @param visualize Logical. If TRUE, saves a TIFF + Rdata of the weight plot.
#' @param group_name Optional character used in plot titles/filenames.
#' @param save_dir Directory to save plots (default "./SoftHybrid_DIA").
#'
#' @return A matrix with soft-hybrid imputed values (same dim as inputs).
#' @export
#' @importFrom ggplot2 %+replace%
#'
#' @examples
#' \dontrun{
#' imp <- soft_hybrid_impute(df_rf, df_minProb, df_raw, visualize = FALSE)
#' }
soft_hybrid_impute <- function(df_rf, df_minProb, df_raw,
                               r0 = NULL, x0 = NULL, a = 10, b = 5, lambda = 0.5,
                               visualize = TRUE, group_name = NULL,
                               save_dir = "./SoftHybrid_DIA") {
  stopifnot(all(dim(df_rf) == dim(df_minProb)))
  stopifnot(all(rownames(df_rf) == rownames(df_minProb)))

  missing_rate   <- rowMeans(is.na(df_raw))
  mean_intensity <- rowMeans(df_raw, na.rm = TRUE)
  miss_cor <- data.frame(missing_rate = missing_rate, mean_abundance = mean_intensity)

  q3 <- stats::quantile(miss_cor$missing_rate, 0.75, na.rm = TRUE)
  skip_loess <- q3 == 0

  if (is.null(r0) || is.null(x0)) {
    if (skip_loess) {
      elbow_x <- 1 / ncol(df_raw)
      candidates <- mean_intensity[round(missing_rate, 6) == round(elbow_x, 6)]
      elbow_y <- if (length(candidates) == 0) {
        mean(mean_intensity, na.rm = TRUE)
      } else {
        mean(candidates, na.rm = TRUE)
      }
    } else {
      loess_fit <- stats::loess(mean_abundance ~ missing_rate, data = miss_cor, span = 0.75)
      x_vals <- seq(min(miss_cor$missing_rate), max(miss_cor$missing_rate), length.out = 1000)
      y_vals <- stats::predict(loess_fit, newdata = data.frame(missing_rate = x_vals))
      x_start <- x_vals[1]; y_start <- y_vals[1]
      x_end   <- x_vals[length(x_vals)]; y_end <- y_vals[length(y_vals)]
      line_vec <- c(x_end - x_start, y_end - y_start)
      distances <- vapply(seq_along(x_vals), function(i) {
        point_vec <- c(x_vals[i] - x_start, y_vals[i] - y_start)
        cross_prod <- abs(line_vec[1] * point_vec[2] - line_vec[2] * point_vec[1])
        norm_line <- sqrt(sum(line_vec^2))
        cross_prod / norm_line
      }, numeric(1))
      max_index <- which.max(distances)
      elbow_x <- x_vals[max_index]
      elbow_y <- y_vals[max_index]
    }
    r0 <- if (is.null(r0)) elbow_x else r0
    x0 <- if (is.null(x0)) elbow_y else x0
  }

  sigmoid <- function(x, k, x0) 1 / (1 + exp(-k * (x - x0)))
  sigmoid_r <- sigmoid(missing_rate, k = a, x0 = r0)
  sigmoid_x <- sigmoid(mean_intensity, k = b, x0 = x0)
  p_MNAR <- sigmoid_r * (1 - lambda * sigmoid_x)

  df_imp_soft <- as.matrix(df_rf)
  for (i in seq_len(nrow(df_rf))) {
    w <- 1 - p_MNAR[i]
    df_imp_soft[i, ] <- w * df_rf[i, ] + (1 - w) * df_minProb[i, ]
  }

  if (isTRUE(visualize)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    .plot_soft_weights(missing_rate, mean_intensity, 1 - p_MNAR, group_name, save_dir)
  }

  df_imp_soft
}

#' Run Soft Hybrid Across Groups
#'
#' @param raw_list A named list of raw matrices.
#' @param rf_list A named list of lists with element \code{imputed_data}.
#' @param minProb_list A named list of lists with element \code{imputed_data}.
#' @inheritParams soft_hybrid_impute
#' @return A named list of matrices with soft-hybrid imputations.
#' @export
run_soft_hybrid_all_groups <- function(raw_list, rf_list, minProb_list,
                                       lambda = 0.5, a = 10, b = 5, visualize = TRUE) {
  group_names <- names(raw_list)
  result_list <- vector("list", length(group_names))
  names(result_list) <- group_names

  for (group in group_names) {
    message("Processing group: ", group)
    df_raw     <- raw_list[[group]]
    df_rf      <- rf_list[[group]][['imputed_data']]
    df_minProb <- minProb_list[[group]][['imputed_data']]

    result_list[[group]] <- soft_hybrid_impute(
      df_rf, df_minProb, df_raw,
      a = a, b = b, lambda = lambda, visualize = visualize,
      group_name = group
    )
  }
  result_list
}

# internal plotting helper
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient scale_x_continuous
#'   scale_y_continuous labs theme_classic guides guide_colorbar
#' @importFrom ggplot2 ggsave
#' @importFrom grid unit
#' @importFrom scales percent_format
.plot_soft_weights <- function(missing_rate, mean_intensity, weight_rf, group_name, save_dir) {
  df_plot <- data.frame(
    MissingRate = missing_rate,
    MeanIntensity = mean_intensity,
    Weight_RF = weight_rf
  )

  theme_pub <- function(base_size = 10, base_family = "Arial") {
    ggplot2::theme_classic(base_size = base_size, base_family = base_family) %+replace%
      ggplot2::theme(
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.title   = ggplot2::element_text(size = base_size + 1, face = "bold"),
        axis.text    = ggplot2::element_text(size = base_size),
        axis.ticks.length = grid::unit(2, "mm"),
        plot.title   = ggplot2::element_text(size = base_size + 2, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title    = ggplot2::element_text(size = base_size),
        legend.text     = ggplot2::element_text(size = base_size - 1)
      )
  }

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = MissingRate, y = MeanIntensity, color = Weight_RF)) +
    ggplot2::geom_point(size = 1.6, alpha = 0.8) +
    ggplot2::scale_color_gradient(low = "blue", high = "red", limits = c(0, 1)) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
    ggplot2::labs(
      title = paste0("Soft Hybrid Weights (", group_name, ")"),
      x = "Missing Rate",
      y = "Mean Intensity",
      color = "RF Weight"
    ) +
    theme_pub(base_size = 10, base_family = "Arial") +
    ggplot2::guides(color = ggplot2::guide_colorbar(
      ticks = TRUE,
      barheight = grid::unit(48, "pt"),
      frame.colour = "black"
    ))

  saveRDS(p, file = file.path(save_dir, paste0(group_name, "_SoftHybrid_plot.rds")))
  ggplot2::ggsave(
    filename = file.path(save_dir, paste0(group_name, "_SoftHybrid_plot.tiff")),
    plot = p, device = "tiff", dpi = 300, width = 6, height = 3.5, units = "in"
  )
  invisible(p)
}
