#' Read Seurat Spatial Objects into SFinxSDFA Format
#'
#' Converts a named list of Seurat spatial transcriptomics objects (one per
#' condition/time point) into a unified \code{data.frame} compatible with
#' \code{STFA_fit()} and other SFinxSDFA functions. Supports all major ST
#' platforms: 10x Visium, Visium HD, Xenium, MERFISH, Slide-seq, and STEREO-seq.
#'
#' @param seurat_list A named list of Seurat objects. The \strong{names} of the
#'   list become the \code{Condition} labels in the output (e.g.,
#'   \code{list(Baseline = seu1, Week4 = seu2, Response = seu3)}).
#' @param pathways Character vector of pathway or gene score column names to
#'   extract from each Seurat object's active assay via \code{FetchData()}.
#'   These must exist as features in the Seurat object's active assay.
#' @param patient_col Optional character string naming a metadata column in each
#'   Seurat object that identifies the patient/donor ID. If \code{NULL}
#'   (default), the condition name is used as \code{PatientID}.
#' @param image_slot Character string naming the image slot to extract from each
#'   Seurat object (default \code{"slice1"}). Used for histology overlays in
#'   downstream plotting. Set to \code{NULL} to skip image extraction.
#'
#' @return A \code{data.frame} with the following columns:
#' \describe{
#'   \item{Longitude}{Numeric x-coordinate of each spot/cell.}
#'   \item{Latitude}{Numeric y-coordinate of each spot/cell.}
#'   \item{Condition}{Factor column with condition/time-point labels from
#'     \code{names(seurat_list)}.}
#'   \item{PatientID}{Character column with patient/donor ID.}
#'   \item{...}{One numeric column per element of \code{pathways}.}
#' }
#' The returned object also has an attribute \code{"images"}: a named list of
#' image matrices (one per condition) for use with \code{plot_stfa_static()}.
#'
#' @examples
#' # --- Synthetic example (no Seurat required) ---
#' # Build mock data directly instead of using Seurat objects
#' set.seed(42)
#' mock_data <- make_st_synthetic(
#'   n_spots = 100,
#'   conditions = c("Baseline", "Week4"),
#'   pathways = c("GOBP_NEUTROPHIL_ACTIVATION", "HP_LUPUS_NEPHRITIS")
#' )
#' # In real use with Seurat objects:
#' # result <- read_seurat_st(
#' #   seurat_list = list(Baseline = seu_baseline, Week4 = seu_week4),
#' #   pathways = c("GOBP_NEUTROPHIL_ACTIVATION", "HP_LUPUS_NEPHRITIS")
#' # )
#'
#' @export
read_seurat_st <- function(seurat_list,
                           pathways,
                           patient_col  = NULL,
                           image_slot   = "slice1") {

  # --- Input validation ---
  if (!is.list(seurat_list) || is.null(names(seurat_list))) {
    stop(
      "[SFinxSDFA] 'seurat_list' must be a named list of Seurat objects. ",
      "Names define the Condition labels (e.g., list(Baseline = seu1, Week4 = seu2))."
    )
  }
  if (length(seurat_list) < 2L) {
    stop(
      "[SFinxSDFA] 'seurat_list' must contain at least 2 conditions ",
      "(got ", length(seurat_list), ")."
    )
  }
  if (!is.character(pathways) || length(pathways) == 0L) {
    stop("[SFinxSDFA] 'pathways' must be a non-empty character vector of feature names.")
  }

  # --- Seurat availability check ---
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(
      "[SFinxSDFA] The 'SeuratObject' package is required to use read_seurat_st(). ",
      "Install it with: install.packages('SeuratObject') or ",
      "remotes::install_github('satijalab/seurat-object')"
    )
  }

  condition_names <- names(seurat_list)
  all_rows        <- vector("list", length(seurat_list))
  images_out      <- vector("list", length(seurat_list))
  names(images_out) <- condition_names

  for (i in seq_along(seurat_list)) {
    cond <- condition_names[i]
    obj  <- seurat_list[[i]]

    # --- Extract spatial coordinates ---
    coords <- tryCatch(
      SeuratObject::GetTissueCoordinates(obj),
      error = function(e) {
        stop(
          "[SFinxSDFA] Could not extract tissue coordinates from Seurat object '", cond, "'. ",
          "Error: ", conditionMessage(e)
        )
      }
    )

    # Harmonize coordinate column names to Longitude / Latitude
    coord_cols <- colnames(coords)
    # Common coordinate column name patterns across platforms
    lon_candidates <- c("imagerow", "x", "col", "array_col", "x_centroid", "spatial_1")
    lat_candidates <- c("imagecol", "y", "row", "array_row", "y_centroid", "spatial_2")

    lon_col <- intersect(tolower(coord_cols), lon_candidates)[1]
    lat_col <- intersect(tolower(coord_cols), lat_candidates)[1]

    # Fall back: use first two numeric columns
    if (is.na(lon_col) || is.na(lat_col)) {
      numeric_cols <- coord_cols[sapply(coords, is.numeric)]
      if (length(numeric_cols) < 2L) {
        stop(
          "[SFinxSDFA] Cannot identify coordinate columns in Seurat object '", cond, "'. ",
          "Columns found: ", paste(coord_cols, collapse = ", ")
        )
      }
      lon_col <- numeric_cols[1]
      lat_col <- numeric_cols[2]
      message(
        "[SFinxSDFA] Warning: using columns '", lon_col, "' and '", lat_col,
        "' as Longitude/Latitude for condition '", cond, "'."
      )
    } else {
      # Match back to original case
      lon_col <- coord_cols[tolower(coord_cols) == lon_col][1]
      lat_col <- coord_cols[tolower(coord_cols) == lat_col][1]
    }

    longitude <- as.numeric(coords[[lon_col]])
    latitude  <- as.numeric(coords[[lat_col]])
    n_spots   <- length(longitude)

    # --- Extract pathway scores ---
    scores <- tryCatch(
      SeuratObject::FetchData(obj, vars = pathways),
      error = function(e) {
        stop(
          "[SFinxSDFA] Could not fetch pathway scores from Seurat object '", cond, "'. ",
          "Check that all 'pathways' exist as features in the active assay. ",
          "Error: ", conditionMessage(e)
        )
      }
    )
    # Ensure row alignment with coordinates
    scores <- scores[rownames(coords), , drop = FALSE]

    # --- Patient ID ---
    if (!is.null(patient_col)) {
      if (!patient_col %in% colnames(obj@meta.data)) {
        stop(
          "[SFinxSDFA] 'patient_col' = '", patient_col, "' not found in metadata of ",
          "Seurat object '", cond, "'. Available columns: ",
          paste(colnames(obj@meta.data), collapse = ", ")
        )
      }
      patient_id <- as.character(obj@meta.data[[patient_col]])
    } else {
      patient_id <- rep(cond, n_spots)
    }

    # --- Assemble per-condition data frame ---
    df_cond <- data.frame(
      Longitude = longitude,
      Latitude  = latitude,
      Condition = cond,
      PatientID = patient_id,
      stringsAsFactors = FALSE
    )
    df_cond <- cbind(df_cond, as.data.frame(scores))

    all_rows[[i]] <- df_cond

    # --- Extract histology image (optional) ---
    if (!is.null(image_slot)) {
      img <- tryCatch({
        img_names <- SeuratObject::Images(obj)
        if (image_slot %in% img_names) {
          SeuratObject::GetImage(obj, image = image_slot, mode = "raster")
        } else if (length(img_names) > 0L) {
          message(
            "[SFinxSDFA] Image slot '", image_slot, "' not found in '", cond,
            "'. Using first available: '", img_names[1], "'."
          )
          SeuratObject::GetImage(obj, image = img_names[1], mode = "raster")
        } else {
          NULL
        }
      }, error = function(e) NULL)

      images_out[[cond]] <- if (is.null(img)) {
        matrix("#FFFFFF", nrow = 1L, ncol = 1L)
      } else {
        img
      }
    } else {
      images_out[[cond]] <- matrix("#FFFFFF", nrow = 1L, ncol = 1L)
    }
  }

  # --- Combine all conditions ---
  result            <- do.call(rbind, all_rows)
  result$Condition  <- factor(result$Condition, levels = condition_names)
  result$PatientID  <- as.character(result$PatientID)
  rownames(result)  <- NULL

  # Attach images as attribute for downstream plotting
  attr(result, "images") <- images_out

  message(
    "[SFinxSDFA] read_seurat_st(): ", nrow(result), " spots across ",
    length(condition_names), " conditions (",
    paste(condition_names, collapse = ", "), ")."
  )

  return(result)
}


#' Create Synthetic Spatiotemporal Data (No Seurat Required)
#'
#' Helper function that generates a synthetic multi-condition
#' \code{data.frame} in SFinxSDFA format for testing and vignettes.
#' Produces biologically realistic spatial patterns with condition-dependent
#' pathway gradients.
#'
#' @param n_spots Integer. Number of spots per condition (default 150).
#' @param conditions Character vector of condition labels
#'   (default \code{c("Baseline", "Week4", "Response")}).
#' @param pathways Character vector of pathway names. Scores are drawn from
#'   multivariate normal distributions with spatially correlated structure.
#' @param seed Integer random seed for reproducibility (default 42).
#'
#' @return A \code{data.frame} with columns \code{Longitude}, \code{Latitude},
#'   \code{Condition} (factor), \code{PatientID}, and one column per pathway.
#'   Has an \code{"images"} attribute (list of white placeholder matrices).
#'
#' @examples
#' df <- make_st_synthetic(n_spots = 100, conditions = c("Pre", "Post"))
#' head(df)
#' table(df$Condition)
#'
#' @export
make_st_synthetic <- function(n_spots   = 150L,
                               conditions = c("Baseline", "Week4", "Response"),
                               pathways   = c(
                                 "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
                                 "HP_LUPUS_NEPHRITIS"
                               ),
                               seed = 42L) {

  set.seed(seed)
  n_cond <- length(conditions)

  all_rows <- vector("list", n_cond)

  for (i in seq_len(n_cond)) {
    cond <- conditions[i]

    lon <- stats::runif(n_spots, 0, 100)
    lat <- stats::runif(n_spots, 0, 100)

    # Introduce a spatial gradient that evolves across conditions
    # so that the spatiotemporal model has real signal to fit
    gradient <- sin(lon / 100 * pi) * cos(lat / 100 * pi) * i

    scores <- lapply(seq_along(pathways), function(j) {
      base  <- stats::rnorm(n_spots, mean = 0, sd = 1.2)
      base + gradient * (0.5 + 0.3 * j)
    })
    scores_df <- as.data.frame(stats::setNames(scores, pathways))

    all_rows[[i]] <- data.frame(
      Longitude = lon,
      Latitude  = lat,
      Condition = cond,
      PatientID = paste0("Patient_", i),
      stringsAsFactors = FALSE
    )
    all_rows[[i]] <- cbind(all_rows[[i]], scores_df)
  }

  result           <- do.call(rbind, all_rows)
  result$Condition <- factor(result$Condition, levels = conditions)
  rownames(result) <- NULL

  # Attach blank image placeholders
  images_out <- stats::setNames(
    lapply(conditions, function(x) matrix("#FFFFFF", nrow = 1L, ncol = 1L)),
    conditions
  )
  attr(result, "images") <- images_out

  return(result)
}
