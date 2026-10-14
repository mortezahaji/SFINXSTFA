#' Spatiotemporal Functional Data Analysis (STFA) — Multi-Condition GAM
#'
#' Fits a spatiotemporally varying coefficient Generalized Additive Model (GAM)
#' across multiple conditions or time points using \code{mgcv::bam} with tensor
#' product P-splines for the spatial surface and Gaussian Process (GP) smooths
#' for condition-specific spatial variation. An optional patient-level random
#' effect accounts for repeated measurements from the same donor.
#'
#' @details
#' The model formula is:
#' \preformatted{
#'   Response ~ te(Longitude, Latitude, bs='ps', k=k)           # global spatial surface
#'            + te(Longitude, Latitude, by=Predictor, bs='ps', k=k)  # spatially varying coeff
#'            + s(Longitude, Latitude, by=Condition, bs='gp')   # condition-specific GP deviation
#'            + s(PatientID, bs='re')                            # patient random effect
#' }
#' Per-condition beta surfaces are extracted by evaluating the fitted model
#' on a \code{grid_n x grid_n} prediction grid with \code{Predictor = 1}
#' (i.e., the spatially varying coefficient surface for that condition).
#' Pointwise 95% confidence intervals are computed from the standard error of
#' the linear predictor. A significance mask is set to 1 where the CI does not
#' cross zero.
#'
#' @param Data A \code{data.frame} with columns:
#'   \code{Longitude}, \code{Latitude}, \code{Condition} (factor),
#'   \code{PatientID}, and the columns named by \code{Response} and
#'   \code{Predictor}.
#' @param Response Character string. Name of the response variable column
#'   (pathway to be predicted).
#' @param Predictor Character string. Name of the predictor variable column
#'   (pathway whose spatially varying effect is estimated).
#' @param k Integer. Basis dimension for tensor product P-splines (default 12).
#'   Increase for smoother surfaces at the cost of computation time.
#' @param grid_n Integer. Resolution of the 2D prediction grid (default 30).
#'   A \code{grid_n x grid_n} grid is evaluated for each condition.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{model}{The fitted \code{bam} model object.}
#'   \item{conditions}{Character vector of condition levels.}
#'   \item{longitude}{Numeric vector of grid longitude values (\code{grid_n} elements).}
#'   \item{latitude}{Numeric vector of grid latitude values (\code{grid_n} elements).}
#'   \item{beta_surfaces}{Named list of \code{grid_n x grid_n} matrices of
#'     spatially varying coefficient estimates per condition.}
#'   \item{lower_ci}{Named list of \code{grid_n x grid_n} lower 95% CI matrices.}
#'   \item{upper_ci}{Named list of \code{grid_n x grid_n} upper 95% CI matrices.}
#'   \item{significance_masks}{Named list of binary \code{grid_n x grid_n} matrices
#'     (1 = significant at 95% level, 0 = not significant).}
#' }
#'
#' @examples
#' # Generate synthetic 3-condition dataset
#' df <- make_st_synthetic(n_spots = 100, conditions = c("Baseline", "Week4", "Response"))
#'
#' # Fit the spatiotemporal model
#' fit <- STFA_fit(
#'   Data      = df,
#'   Response  = "HP_LUPUS_NEPHRITIS",
#'   Predictor = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
#'   k         = 8,
#'   grid_n    = 20
#' )
#'
#' # Inspect outputs
#' names(fit$beta_surfaces)          # condition names
#' dim(fit$beta_surfaces[["Baseline"]])  # grid_n x grid_n
#' summary(fit$model)
#'
#' @export
STFA_fit <- function(Data, Response, Predictor, k = 12L, grid_n = 30L) {

  # --- Input validation ---
  required_cols <- c("Longitude", "Latitude", "Condition", "PatientID",
                     Response, Predictor)
  missing_cols  <- setdiff(required_cols, colnames(Data))
  if (length(missing_cols) > 0L) {
    stop(
      "[SFinxSDFA] STFA_fit(): missing required columns in 'Data': ",
      paste(missing_cols, collapse = ", "), ".\n",
      "  Available columns: ", paste(colnames(Data), collapse = ", ")
    )
  }

  if (!is.factor(Data$Condition)) {
    stop(
      "[SFinxSDFA] STFA_fit(): 'Condition' column must be a factor. ",
      "Convert it with: Data$Condition <- factor(Data$Condition)"
    )
  }

  conditions <- levels(Data$Condition)
  if (length(conditions) < 2L) {
    stop(
      "[SFinxSDFA] STFA_fit(): at least 2 condition levels are required ",
      "(got ", length(conditions), ")."
    )
  }

  n_total <- nrow(Data)
  message(
    "[SFinxSDFA] STFA_fit(): fitting spatiotemporal GAM on ", n_total,
    " spots across ", length(conditions), " conditions (",
    paste(conditions, collapse = ", "), ")."
  )
  message(
    "[SFinxSDFA]   Response: '", Response, "'  |  Predictor: '", Predictor,
    "'  |  k = ", k, "  |  grid_n = ", grid_n
  )

  # --- Check for patient random effect viability ---
  # If all PatientID values are unique (one patient per condition) we can still
  # fit the model but the re smooth will degenerate; warn but proceed.
  n_patients <- length(unique(Data$PatientID))
  use_re     <- n_patients > 1L
  if (!use_re) {
    message(
      "[SFinxSDFA] Warning: only one unique PatientID found. ",
      "Patient random effect (bs='re') will be omitted."
    )
  }

  # --- Build GAM formula ---
  # Each condition level gets its own GP smooth for spatial deviation.
  # We create a separate s() term per condition level by expanding them
  # explicitly (mgcv requires factor 'by' variable to be ordered correctly).
  cond_smooth_terms <- paste0(
    "s(Longitude, Latitude, bs='gp', by=",
    # Use interaction of Condition dummy with the predictor
    # We'll dummy-encode conditions inline via the formula
    "I(as.numeric(Condition == '", conditions, "')))",
    collapse = " + "
  )

  re_term <- if (use_re) " + s(PatientID, bs='re')" else ""

  formula_str <- paste0(
    Response, " ~ ",
    "te(Longitude, Latitude, bs='ps', k=", k, ") + ",
    "te(Longitude, Latitude, by=", Predictor, ", bs='ps', k=", k, ") + ",
    cond_smooth_terms,
    re_term
  )

  message("[SFinxSDFA] Model formula: ", formula_str)

  # Make PatientID a factor for the random effect
  Data$PatientID <- factor(Data$PatientID)

  # --- Fit model ---
  fit_model <- tryCatch(
    mgcv::bam(
      stats::as.formula(formula_str),
      data   = Data,
      method = "REML"
    ),
    error = function(e) {
      stop(
        "[SFinxSDFA] STFA_fit(): model fitting failed.\n",
        "  Error: ", conditionMessage(e), "\n",
        "  Try reducing 'k' or increasing 'n_spots' per condition."
      )
    }
  )

  message("[SFinxSDFA] Model fitting complete. Extracting per-condition surfaces...")

  # --- Build prediction grid ---
  lon_vals  <- seq(min(Data$Longitude), max(Data$Longitude), length.out = grid_n)
  lat_vals  <- seq(min(Data$Latitude),  max(Data$Latitude),  length.out = grid_n)
  base_grid <- expand.grid(Longitude = lon_vals, Latitude = lat_vals)

  # --- Extract per-condition beta surfaces ---
  beta_surfaces     <- vector("list", length(conditions))
  lower_ci_list     <- vector("list", length(conditions))
  upper_ci_list     <- vector("list", length(conditions))
  sig_mask_list     <- vector("list", length(conditions))
  names(beta_surfaces) <- conditions
  names(lower_ci_list) <- conditions
  names(upper_ci_list) <- conditions
  names(sig_mask_list) <- conditions

  # Reference PatientID for prediction (use first level)
  ref_patient <- levels(Data$PatientID)[1]

  for (cond in conditions) {
    pred_grid <- base_grid

    # Set Predictor = 1 to extract the spatially varying coefficient surface
    pred_grid[[Predictor]]  <- 1
    pred_grid[[Response]]   <- 0   # placeholder (not used in prediction)

    # Encode condition dummies
    for (c2 in conditions) {
      dummy_col <- paste0("I(as.numeric(Condition == '", c2, "'))")
      pred_grid[[dummy_col]] <- as.numeric(c2 == cond)
    }

    # Set Condition and PatientID for prediction
    pred_grid$Condition <- factor(cond, levels = conditions)
    pred_grid$PatientID <- factor(ref_patient, levels = levels(Data$PatientID))

    # Predict with SE
    pred_obj <- tryCatch(
      stats::predict(fit_model, newdata = pred_grid, type = "link", se.fit = TRUE),
      error = function(e) {
        message(
          "[SFinxSDFA] Warning: prediction failed for condition '", cond,
          "'. Returning zero surface. Error: ", conditionMessage(e)
        )
        list(
          fit   = rep(0, grid_n * grid_n),
          se.fit = rep(0, grid_n * grid_n)
        )
      }
    )

    beta_vec  <- as.numeric(pred_obj$fit)
    se_vec    <- as.numeric(pred_obj$se.fit)
    lower_vec <- beta_vec - 1.96 * se_vec
    upper_vec <- beta_vec + 1.96 * se_vec
    sig_vec   <- as.numeric((lower_vec > 0) | (upper_vec < 0))

    beta_surfaces[[cond]] <- matrix(beta_vec,  nrow = grid_n, ncol = grid_n)
    lower_ci_list[[cond]] <- matrix(lower_vec, nrow = grid_n, ncol = grid_n)
    upper_ci_list[[cond]] <- matrix(upper_vec, nrow = grid_n, ncol = grid_n)
    sig_mask_list[[cond]] <- matrix(sig_vec,   nrow = grid_n, ncol = grid_n)

    message("[SFinxSDFA]   Condition '", cond, "' surface extracted.")
  }

  message("[SFinxSDFA] STFA_fit() complete.")

  return(list(
    model              = fit_model,
    conditions         = conditions,
    longitude          = lon_vals,
    latitude           = lat_vals,
    beta_surfaces      = beta_surfaces,
    lower_ci           = lower_ci_list,
    upper_ci           = upper_ci_list,
    significance_masks = sig_mask_list
  ))
}
