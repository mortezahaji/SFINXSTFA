#' Interactive Plotly Spatial Explorer with Condition Slider
#'
#' Produces an interactive HTML figure from \code{\link{STFA_fit}} output.
#' Pathologists can scrub through conditions using a slider and toggle between
#' a 2D heatmap view (with significance contour overlay) and a 3D surface view
#' of the spatially varying beta coefficient surface.
#'
#' @details
#' The returned \code{plotly} object contains one trace per condition, with all
#' but the first initially hidden. A layout slider controls visibility. The
#' slider labels show the condition names so pathologists see e.g. "Baseline",
#' "Week 4", "Response" rather than index numbers.
#'
#' @param stfa_result A list returned by \code{\link{STFA_fit}}.
#' @param title Character string. Title displayed above the figure
#'   (default \code{"Spatiotemporal Pathway Explorer"}).
#' @param mode Character string. Either \code{"2d"} (default) for a 2D heatmap
#'   with significance contour overlay, or \code{"3d"} for a 3D surface.
#'
#' @return A \code{plotly} / \code{htmlwidget} object. Save with:
#' \preformatted{
#'   htmlwidgets::saveWidget(p, "explorer.html")
#'   # or
#'   plotly::save_image(p, "explorer.png")
#' }
#'
#' @examples
#' df  <- make_st_synthetic(n_spots = 80, conditions = c("Baseline", "Week4", "Response"))
#' fit <- STFA_fit(df, Response = "HP_LUPUS_NEPHRITIS",
#'                 Predictor = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
#'                 k = 6, grid_n = 15)
#'
#' # 2D heatmap with slider
#' p2d <- plot_stfa_interactive(fit, mode = "2d")
#' # print(p2d)  # opens in browser / RStudio viewer
#'
#' # 3D surface with slider
#' p3d <- plot_stfa_interactive(fit, mode = "3d")
#' # print(p3d)
#'
#' @export
plot_stfa_interactive <- function(stfa_result,
                                   title = "Spatiotemporal Pathway Explorer",
                                   mode  = c("2d", "3d")) {

  mode <- match.arg(mode)

  # --- Validate input ---
  required_elements <- c("conditions", "longitude", "latitude",
                         "beta_surfaces", "significance_masks")
  missing_el <- setdiff(required_elements, names(stfa_result))
  if (length(missing_el) > 0L) {
    stop(
      "[SFinxSDFA] plot_stfa_interactive(): 'stfa_result' is missing elements: ",
      paste(missing_el, collapse = ", "),
      ". Pass the output of STFA_fit() directly."
    )
  }

  conditions <- stfa_result$conditions
  lon_vals   <- stfa_result$longitude
  lat_vals   <- stfa_result$latitude
  n_cond     <- length(conditions)

  # Global beta range for shared colour axis
  all_betas  <- unlist(lapply(stfa_result$beta_surfaces, as.numeric))
  beta_range <- range(all_betas, na.rm = TRUE)

  # --- Build plotly figure ---
  if (mode == "2d") {
    fig <- .build_2d_slider(
      stfa_result = stfa_result,
      conditions  = conditions,
      lon_vals    = lon_vals,
      lat_vals    = lat_vals,
      beta_range  = beta_range,
      title       = title
    )
  } else {
    fig <- .build_3d_slider(
      stfa_result = stfa_result,
      conditions  = conditions,
      lon_vals    = lon_vals,
      lat_vals    = lat_vals,
      beta_range  = beta_range,
      title       = title
    )
  }

  return(fig)
}


# ---------------------------------------------------------------------------
# Internal helper: 2D heatmap with significance contour overlay + slider
# ---------------------------------------------------------------------------
.build_2d_slider <- function(stfa_result, conditions, lon_vals, lat_vals,
                              beta_range, title) {

  n_cond        <- length(conditions)
  traces_per    <- 2L  # heatmap + significance contour per condition
  total_traces  <- n_cond * traces_per

  # Visibility matrix: each slider step shows only its own pair of traces
  vis_matrix <- lapply(seq_len(n_cond), function(active_i) {
    vis <- rep(FALSE, total_traces)
    vis[((active_i - 1L) * traces_per + 1L):(active_i * traces_per)] <- TRUE
    vis
  })

  fig <- plotly::plot_ly()

  for (i in seq_len(n_cond)) {
    cond    <- conditions[i]
    beta_m  <- stfa_result$beta_surfaces[[cond]]
    sig_m   <- stfa_result$significance_masks[[cond]]
    visible <- (i == 1L)

    # --- Heatmap trace ---
    fig <- plotly::add_trace(
      fig,
      x          = lon_vals,
      y          = lat_vals,
      z          = beta_m,
      type       = "heatmap",
      colorscale = "Viridis",
      zmin       = beta_range[1],
      zmax       = beta_range[2],
      colorbar   = list(
        title      = "\u03b2(s)",
        titlefont  = list(size = 13),
        tickfont   = list(size = 11)
      ),
      name       = paste0(cond, " \u03b2"),
      visible    = visible,
      showscale  = visible,
      hovertemplate = paste0(
        "<b>", cond, "</b><br>",
        "Longitude: %{x:.1f}<br>",
        "Latitude: %{y:.1f}<br>",
        "\u03b2: %{z:.3f}<extra></extra>"
      )
    )

    # --- Significance contour overlay trace ---
    # Create a contour plot of the significance mask (shows where CI ≠ 0)
    fig <- plotly::add_trace(
      fig,
      x          = lon_vals,
      y          = lat_vals,
      z          = sig_m,
      type       = "contour",
      contours   = list(
        start = 0.5,
        end   = 1.5,
        size  = 1.0,
        showlabels = FALSE
      ),
      line       = list(color = "white", width = 1.5, dash = "dot"),
      showscale  = FALSE,
      name       = paste0(cond, " sig"),
      visible    = visible,
      hoverinfo  = "skip"
    )
  }

  # --- Build slider steps ---
  steps <- lapply(seq_len(n_cond), function(i) {
    list(
      args   = list("visible", vis_matrix[[i]]),
      label  = conditions[i],
      method = "restyle"
    )
  })

  fig <- plotly::layout(
    fig,
    title = list(
      text = title,
      font = list(size = 15, family = "Arial")
    ),
    xaxis  = list(title = "Longitude", zeroline = FALSE),
    yaxis  = list(title = "Latitude",  zeroline = FALSE),
    sliders = list(list(
      active      = 0L,
      currentvalue = list(
        prefix = "Condition: ",
        font   = list(size = 12, color = "#444")
      ),
      steps = steps,
      pad   = list(t = 50L)
    )),
    margin = list(t = 100L, b = 80L),
    annotations = list(list(
      text      = paste0(
        "White dotted contour = significant (95% CI \u2260 0)"
      ),
      xref      = "paper",
      yref      = "paper",
      x         = 0,
      y         = -0.12,
      showarrow = FALSE,
      font      = list(size = 10, color = "grey50")
    ))
  )

  return(fig)
}


# ---------------------------------------------------------------------------
# Internal helper: 3D surface with slider
# ---------------------------------------------------------------------------
.build_3d_slider <- function(stfa_result, conditions, lon_vals, lat_vals,
                              beta_range, title) {

  n_cond <- length(conditions)

  # For 3D, each condition = 1 surface trace
  # Visibility: one surface visible at a time
  vis_matrix <- lapply(seq_len(n_cond), function(active_i) {
    vis <- rep(FALSE, n_cond)
    vis[active_i] <- TRUE
    vis
  })

  fig <- plotly::plot_ly()

  for (i in seq_len(n_cond)) {
    cond    <- conditions[i]
    beta_m  <- stfa_result$beta_surfaces[[cond]]
    visible <- (i == 1L)

    fig <- plotly::add_surface(
      fig,
      x          = lon_vals,
      y          = lat_vals,
      z          = beta_m,
      colorscale = "Viridis",
      cmin       = beta_range[1],
      cmax       = beta_range[2],
      contours   = list(
        z = list(
          show         = TRUE,
          usecolormap  = TRUE,
          highlightcolor = "#ffffff",
          project      = list(z = TRUE)
        )
      ),
      colorbar  = list(
        title    = "\u03b2(s)",
        titlefont = list(size = 13)
      ),
      name      = cond,
      visible   = visible,
      showscale = visible
    )
  }

  # --- Slider steps ---
  steps <- lapply(seq_len(n_cond), function(i) {
    list(
      args   = list("visible", vis_matrix[[i]]),
      label  = conditions[i],
      method = "restyle"
    )
  })

  fig <- plotly::layout(
    fig,
    title = list(
      text = title,
      font = list(size = 15, family = "Arial")
    ),
    scene = list(
      xaxis = list(title = "Longitude"),
      yaxis = list(title = "Latitude"),
      zaxis = list(title = "\u03b2(s)")
    ),
    sliders = list(list(
      active       = 0L,
      currentvalue = list(
        prefix = "Condition: ",
        font   = list(size = 12, color = "#444")
      ),
      steps = steps,
      pad   = list(t = 50L)
    )),
    margin = list(t = 100L, b = 80L)
  )

  return(fig)
}
