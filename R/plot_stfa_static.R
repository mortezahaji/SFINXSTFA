#' Static Faceted Spatial Maps of Spatiotemporal Pathway Analysis
#'
#' Generates a publication-ready multi-panel \code{ggplot2} figure from the
#' output of \code{\link{STFA_fit}}. Each panel shows the spatially varying
#' coefficient (beta) surface for one condition, with a significance overlay
#' and optional histology background. All panels share the same colour scale
#' and axis limits for direct visual comparison across conditions.
#'
#' @param stfa_result A list returned by \code{\link{STFA_fit}}.
#' @param title Character string. Main title displayed above the figure
#'   (default \code{"Spatiotemporal Pathway Analysis"}).
#' @param images Optional named list of image matrices, one per condition
#'   (names must match \code{stfa_result$conditions}). Typically the
#'   \code{attr(data, "images")} attribute from \code{\link{read_seurat_st}} or
#'   \code{\link{make_st_synthetic}}. If \code{NULL}, a plain white background
#'   is used (consistent with \code{SFinx_UO}/\code{SFinx_BO}).
#' @param sig_overlay Logical. Whether to draw the significance overlay as
#'   white semi-transparent points on grid cells where the 95\% CI excludes
#'   zero (default \code{TRUE}).
#' @param sig_alpha Numeric in [0,1]. Alpha transparency of significance
#'   overlay points (default \code{0.35}).
#' @param ncol Integer. Number of columns in the facet grid. Defaults to the
#'   number of conditions (all in one row). Set to e.g. \code{2} for layouts
#'   with many conditions.
#'
#' @return A \code{ggplot} object. Export with
#'   \code{ggplot2::ggsave("output.png", plot = p, width = 14, height = 5, dpi = 300)}.
#'
#' @examples
#' # Fit model on synthetic data
#' df  <- make_st_synthetic(n_spots = 80, conditions = c("Baseline", "Week4", "Response"))
#' fit <- STFA_fit(df, Response = "HP_LUPUS_NEPHRITIS",
#'                 Predictor = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
#'                 k = 6, grid_n = 15)
#'
#' # Create static plot
#' p <- plot_stfa_static(fit, title = "Lupus Nephritis Pathway Dynamics")
#' print(p)
#'
#' # Export
#' # ggplot2::ggsave("ST_analysis.png", plot = p, width = 14, height = 5, dpi = 300)
#'
#' @export
plot_stfa_static <- function(stfa_result,
                              title       = "Spatiotemporal Pathway Analysis",
                              images      = NULL,
                              sig_overlay = TRUE,
                              sig_alpha   = 0.35,
                              ncol        = NULL) {

  # --- Validate input ---
  required_elements <- c("conditions", "longitude", "latitude",
                         "beta_surfaces", "significance_masks")
  missing_el <- setdiff(required_elements, names(stfa_result))
  if (length(missing_el) > 0L) {
    stop(
      "[SFinxSDFA] plot_stfa_static(): 'stfa_result' is missing elements: ",
      paste(missing_el, collapse = ", "),
      ". Pass the output of STFA_fit() directly."
    )
  }

  conditions <- stfa_result$conditions
  lon_vals   <- stfa_result$longitude
  lat_vals   <- stfa_result$latitude
  n_cond     <- length(conditions)
  ncol       <- if (is.null(ncol)) n_cond else as.integer(ncol)

  # Fallback images
  if (is.null(images)) {
    images <- stats::setNames(
      lapply(conditions, function(x) matrix("#FFFFFF", nrow = 1L, ncol = 1L)),
      conditions
    )
  }

  # --- Build long data frame for all conditions ---
  grid_df_list <- vector("list", n_cond)

  for (cond in conditions) {
    beta_mat <- stfa_result$beta_surfaces[[cond]]
    sig_mat  <- stfa_result$significance_masks[[cond]]
    grid_n   <- nrow(beta_mat)

    base_grid <- expand.grid(Longitude = lon_vals, Latitude = lat_vals)
    base_grid$Beta        <- as.numeric(beta_mat)
    base_grid$Significant <- as.integer(sig_mat)
    base_grid$Condition   <- cond

    grid_df_list[[match(cond, conditions)]] <- base_grid
  }

  plot_df <- do.call(rbind, grid_df_list)
  plot_df$Condition <- factor(plot_df$Condition, levels = conditions)

  # Shared beta range for consistent colour scale across facets
  beta_range <- range(plot_df$Beta, na.rm = TRUE)

  # Significance overlay points (only significant grid cells)
  if (sig_overlay) {
    sig_df <- plot_df[plot_df$Significant == 1L, ]
  }

  # --- Shared theme (mirrors SFinx_UO / SFinx_BO) ---
  base_theme <- ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, size = 13, face = "bold"),
      strip.text       = ggplot2::element_text(size = 11, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      axis.title       = ggplot2::element_text(size = 10),
      axis.text        = ggplot2::element_text(size = 9),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(size = 10),
      panel.spacing    = ggplot2::unit(0.8, "lines")
    )

  # --- Build plot ---
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Longitude, y = Latitude, fill = Beta)
  ) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_viridis_c(
      option = "plasma",
      name   = expression(beta(s)),
      limits = beta_range,
      oob    = scales::squish
    ) +
    ggplot2::facet_wrap(~ Condition, ncol = ncol) +
    ggplot2::labs(
      title = title,
      x     = "Longitude",
      y     = "Latitude"
    ) +
    base_theme

  # Significance overlay
  if (sig_overlay && nrow(sig_df) > 0L) {
    p <- p +
      ggplot2::geom_point(
        data    = sig_df,
        mapping = ggplot2::aes(x = Longitude, y = Latitude),
        inherit.aes = FALSE,
        colour  = "white",
        size    = 0.4,
        alpha   = sig_alpha,
        shape   = 16L
      )
  }

  return(p)
}
