#' Run the Full SFinxSDFA Spatiotemporal Pipeline
#'
#' Single entry-point function for pathologists. Accepts either a named list
#' of Seurat objects or a pre-formatted \code{data.frame} and runs the
#' complete spatiotemporal analysis pipeline:
#'
#' \enumerate{
#'   \item \strong{Data ingestion:} \code{\link{read_seurat_st}} (if Seurat
#'     objects provided) or direct use of the supplied \code{data} frame.
#'   \item \strong{Per-condition spatial fingerprints:} Calls
#'     \code{\link[SFinxSDFA]{SFinx_UO}} on each condition to validate the
#'     spatial signal before fitting.
#'   \item \strong{Spatiotemporal GAM:} \code{\link{STFA_fit}} fits a
#'     multi-condition spatially varying coefficient model.
#'   \item \strong{Visualisation:} \code{\link{plot_stfa_static}} produces a
#'     publication-ready faceted map; \code{\link{plot_stfa_interactive}}
#'     produces an interactive HTML explorer.
#'   \item \strong{Optional export:} Static PNG and interactive HTML are saved
#'     to \code{output_dir} if provided.
#' }
#'
#' @param seurat_list A named list of Seurat objects (one per condition), e.g.
#'   \code{list(Baseline = seu1, Week4 = seu2)}. Provide either this argument
#'   or \code{data}, not both.
#' @param data A \code{data.frame} already in SFinxSDFA format (with columns
#'   \code{Longitude}, \code{Latitude}, \code{Condition} (factor),
#'   \code{PatientID}, and pathway columns). Use this if you have already run
#'   \code{\link{read_seurat_st}} or \code{\link{make_st_synthetic}}.
#' @param response Character string. Name of the response pathway column.
#' @param predictor Character string. Name of the predictor pathway column.
#' @param pathways Character vector of pathway names to extract from
#'   \code{seurat_list} (only used when \code{seurat_list} is provided).
#'   Must contain both \code{response} and \code{predictor}.
#' @param k Integer. Spline basis dimension for \code{\link{STFA_fit}}
#'   (default 12).
#' @param grid_n Integer. Prediction grid resolution for
#'   \code{\link{STFA_fit}} (default 30).
#' @param output_dir Optional character string. Path to a directory where
#'   outputs will be saved:
#'   \itemize{
#'     \item \code{<title>_static.png} — faceted static map (300 dpi).
#'     \item \code{<title>_interactive_2d.html} — 2D slider explorer.
#'     \item \code{<title>_interactive_3d.html} — 3D surface explorer.
#'   }
#'   The directory is created if it does not exist.
#' @param title Character string. Title used for plots and output filenames
#'   (default \code{"Spatiotemporal_Analysis"}). Spaces are replaced with
#'   underscores in filenames.
#' @param ... Additional arguments passed to \code{\link{read_seurat_st}}
#'   (e.g. \code{patient_col}, \code{image_slot}).
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{data}{The unified \code{data.frame} used for modelling.}
#'   \item{model_fit}{The list returned by \code{\link{STFA_fit}}, including
#'     the fitted model, per-condition beta surfaces, CIs, and significance masks.}
#'   \item{static_plot}{A \code{ggplot} object (faceted static map).}
#'   \item{interactive_plot}{A \code{plotly} object (2D heatmap slider).}
#' }
#'
#' @examples
#' # --- Synthetic data example (no Seurat required) ---
#' df <- make_st_synthetic(
#'   n_spots    = 150,
#'   conditions = c("Baseline", "Week4", "Response"),
#'   pathways   = c("GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION", "HP_LUPUS_NEPHRITIS")
#' )
#'
#' result <- run_spatiotemporal(
#'   data      = df,
#'   response  = "HP_LUPUS_NEPHRITIS",
#'   predictor = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
#'   k         = 8,
#'   grid_n    = 20,
#'   title     = "Lupus Nephritis Dynamics"
#' )
#'
#' # Inspect outputs
#' print(result$static_plot)
#' # print(result$interactive_plot)  # opens in browser
#' names(result$model_fit$beta_surfaces)
#'
#' @export
run_spatiotemporal <- function(seurat_list = NULL,
                                data        = NULL,
                                response,
                                predictor,
                                pathways    = NULL,
                                k           = 12L,
                                grid_n      = 30L,
                                output_dir  = NULL,
                                title       = "Spatiotemporal_Analysis",
                                ...) {

  # -------------------------------------------------------------------------
  # 1. Input validation
  # -------------------------------------------------------------------------
  if (is.null(seurat_list) && is.null(data)) {
    stop(
      "[SFinxSDFA] run_spatiotemporal(): provide either 'seurat_list' ",
      "(a named list of Seurat objects) or 'data' (a pre-formatted data.frame)."
    )
  }
  if (!is.null(seurat_list) && !is.null(data)) {
    stop(
      "[SFinxSDFA] run_spatiotemporal(): provide only one of 'seurat_list' ",
      "or 'data', not both."
    )
  }
  if (missing(response) || !is.character(response) || length(response) != 1L) {
    stop("[SFinxSDFA] run_spatiotemporal(): 'response' must be a single character string.")
  }
  if (missing(predictor) || !is.character(predictor) || length(predictor) != 1L) {
    stop("[SFinxSDFA] run_spatiotemporal(): 'predictor' must be a single character string.")
  }

  # -------------------------------------------------------------------------
  # 2. Data ingestion
  # -------------------------------------------------------------------------
  if (!is.null(seurat_list)) {
    message("[SFinxSDFA] Reading Seurat data...")

    # Ensure response and predictor are included in pathways
    pathways <- unique(c(pathways, response, predictor))

    combined_data <- read_seurat_st(
      seurat_list = seurat_list,
      pathways    = pathways,
      ...
    )
  } else {
    combined_data <- data
    message(
      "[SFinxSDFA] Using supplied data.frame: ", nrow(combined_data),
      " spots."
    )
  }

  # Validate required columns
  required_cols <- c("Longitude", "Latitude", "Condition", "PatientID",
                     response, predictor)
  missing_cols  <- setdiff(required_cols, colnames(combined_data))
  if (length(missing_cols) > 0L) {
    stop(
      "[SFinxSDFA] run_spatiotemporal(): 'data' is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Ensure Condition is a factor with ≥2 levels
  if (!is.factor(combined_data$Condition)) {
    combined_data$Condition <- factor(combined_data$Condition)
  }
  n_conditions <- nlevels(combined_data$Condition)
  if (n_conditions < 2L) {
    stop(
      "[SFinxSDFA] run_spatiotemporal(): at least 2 condition levels are required ",
      "(got ", n_conditions, ")."
    )
  }

  # -------------------------------------------------------------------------
  # 3. Per-condition spatial fingerprint validation (SFinx_UO)
  # -------------------------------------------------------------------------
  message("[SFinxSDFA] Computing per-condition spatial fingerprints (SFinx_UO)...")

  conditions    <- levels(combined_data$Condition)
  images_attr   <- attr(combined_data, "images")

  for (cond in conditions) {
    cond_df  <- combined_data[combined_data$Condition == cond, ]
    scores   <- cond_df[, c(response, predictor), drop = FALSE]
    pos_df   <- data.frame(
      Longitude = cond_df$Longitude,
      Latitude  = cond_df$Latitude,
      x         = cond_df$Longitude,
      y         = cond_df$Latitude
    )
    img <- if (!is.null(images_attr) && cond %in% names(images_attr)) {
      images_attr[[cond]]
    } else {
      matrix("#FFFFFF", nrow = 1L, ncol = 1L)
    }

    tryCatch({
      # Run SFinx_UO silently to validate spatial signal per condition
      suppressMessages(suppressWarnings(
        SFinx_UO(
          Data   = scores,
          Finpre = response,
          Thresh = 1.0,
          Pos    = pos_df,
          Image  = img,
          Title  = paste0(cond, ": ", response)
        )
      ))
      message("[SFinxSDFA]   Condition '", cond, "': spatial signal validated.")
    }, error = function(e) {
      message(
        "[SFinxSDFA]   Warning: SFinx_UO failed for condition '", cond,
        "' (continuing). Error: ", conditionMessage(e)
      )
    })
  }

  # -------------------------------------------------------------------------
  # 4. Fit spatiotemporal model
  # -------------------------------------------------------------------------
  message("[SFinxSDFA] Fitting spatiotemporal model...")

  model_fit <- STFA_fit(
    Data      = combined_data,
    Response  = response,
    Predictor = predictor,
    k         = k,
    grid_n    = grid_n
  )

  # -------------------------------------------------------------------------
  # 5. Generate visualisations
  # -------------------------------------------------------------------------
  message("[SFinxSDFA] Generating visualizations...")

  images_list <- attr(combined_data, "images")

  static_plot <- plot_stfa_static(
    stfa_result = model_fit,
    title       = gsub("_", " ", title),
    images      = images_list
  )

  interactive_plot <- plot_stfa_interactive(
    stfa_result = model_fit,
    title       = gsub("_", " ", title),
    mode        = "2d"
  )

  interactive_plot_3d <- plot_stfa_interactive(
    stfa_result = model_fit,
    title       = gsub("_", " ", title),
    mode        = "3d"
  )

  # -------------------------------------------------------------------------
  # 6. Optional file export
  # -------------------------------------------------------------------------
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message("[SFinxSDFA] Created output directory: ", output_dir)
    }

    # Sanitise filename
    file_stem <- gsub("[[:space:]]", "_", title)
    file_stem <- gsub("[^A-Za-z0-9_\\-]", "", file_stem)

    # Static PNG
    static_path <- file.path(output_dir, paste0(file_stem, "_static.png"))
    ggplot2::ggsave(
      filename = static_path,
      plot     = static_plot,
      width    = 5 * length(conditions),
      height   = 5,
      dpi      = 300,
      limitsize = FALSE
    )
    message("[SFinxSDFA] Static plot saved: ", static_path)

    # Interactive 2D HTML
    html_2d_path <- file.path(output_dir, paste0(file_stem, "_interactive_2d.html"))
    htmlwidgets::saveWidget(
      widget    = interactive_plot,
      file      = html_2d_path,
      selfcontained = TRUE
    )
    message("[SFinxSDFA] Interactive 2D plot saved: ", html_2d_path)

    # Interactive 3D HTML
    html_3d_path <- file.path(output_dir, paste0(file_stem, "_interactive_3d.html"))
    htmlwidgets::saveWidget(
      widget    = interactive_plot_3d,
      file      = html_3d_path,
      selfcontained = TRUE
    )
    message("[SFinxSDFA] Interactive 3D plot saved: ", html_3d_path)
  }

  message("[SFinxSDFA] Done.")

  return(list(
    data             = combined_data,
    model_fit        = model_fit,
    static_plot      = static_plot,
    interactive_plot = interactive_plot
  ))
}
