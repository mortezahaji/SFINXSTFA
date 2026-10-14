# SFinxSDFA: Spatial Fingerprint Analytics & Spatial Functional Data Analysis

[![R-CMD-check](https://img.shields.io/badge/R--package-SFinxSDFA-blue.svg)](https://github.com/mortezahaji/SFinxSDFA)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.07.15.738836-red.svg)](https://doi.org/10.64898/2026.07.15.738836)

**`SFinxSDFA`** is an R package for reconstructing continuous molecular landscapes,
evaluating spatial pathway operators, and fitting spatially varying regression models on
spatial transcriptomics (ST) data. By combining **Spatial Fingerprint Analytics (SFinx)**
and **Spatial Functional Data Analysis (SDFA)**, this platform moves beyond gene- or
domain-level clustering to model continuous 2D/3D pathway activity and pathway–phenotype
crosstalk across complex tissue sections.

---

## Key Capabilities

- **Spatial Fingerprints (SFinx):** Quantitative, programmable representations of molecular
  signatures across 2D tissue spaces.
- **Spatial Operators:**
  - **Unary Operators:** Spatial autocorrelation (Moran's $I$), activation scaling
    ($|Z| \ge \text{threshold}$), and 2D kernel density estimations.
  - **Binary Operators:** Union, Intersection (Jaccard Index), Difference, and Maximum
    (Dominance) maps between pathway signatures.
- **Spatial Functional Data Analysis (SDFA):** Continuous 2D surface modeling and
  Generalized Additive Models (GAMs) using tensor product B-splines (`mgcv::bam`).
- **Multi-Slice 3D Integration:** Evaluates spatial relationships across consecutive tissue
  sections while accounting for slice-specific random effects along the z-axis.
- **Spatiotemporal Modeling _(v0.2.0)_:** Multi-condition longitudinal analysis tracking
  *where* and *how* pathway activity changes across clinical time points — with Seurat
  support for all major ST platforms.

---

## Installation

```r
# Install remotes if not already available
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install SFinxSDFA
remotes::install_github("mortezahaji/SFinxSDFA")
```

---

## Spatiotemporal Modeling _(v0.2.0)_

The spatiotemporal extension enables pathologists and computational biologists to compare
multi-condition ST data (e.g., Baseline → Week 4 → Treatment Response) through a single
pipeline call. Internally it fits a spatiotemporally varying coefficient GAM using
`mgcv::bam` with tensor product P-splines for the global spatial surface, Gaussian Process
(GP) smooths for condition-specific spatial deviations, and patient-level random effects.

### Quick Start — Synthetic Data (no Seurat required)

```r
library(SFinxSDFA)
library(ggplot2)

# 1. Generate synthetic 3-condition dataset
df <- make_st_synthetic(
  n_spots    = 150,
  conditions = c("Baseline", "Week4", "Response"),
  pathways   = c(
    "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
    "HP_LUPUS_NEPHRITIS"
  )
)

# 2. Run the full pipeline in one call
result <- run_spatiotemporal(
  data      = df,
  response  = "HP_LUPUS_NEPHRITIS",
  predictor = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
  k         = 8,
  grid_n    = 25,
  title     = "Lupus Nephritis Dynamics"
)

# 3. View publication-ready static map
print(result$static_plot)

# 4. Launch interactive HTML explorer
print(result$interactive_plot)   # opens in browser / RStudio viewer
```

### With Real Seurat Objects

```r
# Works with Visium, Xenium, MERFISH, Slide-seq — platform auto-detected
result <- run_spatiotemporal(
  seurat_list = list(
    Baseline = seu_baseline,
    Week4    = seu_week4,
    Response = seu_response
  ),
  response    = "HP_LUPUS_NEPHRITIS",
  predictor   = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
  pathways    = c("HP_LUPUS_NEPHRITIS", "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION"),
  k           = 12,
  grid_n      = 30,
  output_dir  = "./results",
  title       = "Cohort_Analysis"
)
```

### Spatiotemporal API

| Function | Description |
|---|---|
| `make_st_synthetic()` | Generate synthetic multi-condition ST data |
| `read_seurat_st()` | Import Seurat objects → SFinxSDFA data frame |
| `STFA_fit()` | Fit spatiotemporal GAM; returns per-condition β surfaces |
| `plot_stfa_static()` | Faceted ggplot2 condition maps (publication-ready) |
| `plot_stfa_interactive()` | Plotly slider explorer (2D heatmap or 3D surface) |
| `run_spatiotemporal()` | One-call pipeline orchestrator |

See `vignette("spatiotemporal_workflow")` for a full illustrated walkthrough.

---

## Complete End-to-End Example (v0.1 Core Functions)

The following script runs a synthetic spatial transcriptomics workflow covering:

- **Compatibility patching (ggpubr & ggplot2 bindings).**
- **Synthetic coordinate (Longitude, Latitude) and pathway score generation.**
- **SFinx Unary Operator (SFinx_UO).**
- **SFinx Binary Operators (SFinx_BO).**
- **SDFA Spatial Functional Regression (FOF_fit).**
- **3D Interactive & Static Coefficient Surface Plotting.**

### 1. Load Libraries & Apply Patches

```r
library(SFinxSDFA)
library(ggplot2)
library(ggpubr)
library(plotly)
library(mgcv)

set.seed(42)

# Patch hardcoded namespace calls for local session compatibility
body_UO <- deparse(body(SFinx_UO))
body_UO_fixed <- gsub("egg::background_image", "ggpubr::background_image", body_UO)
body_UO_fixed <- gsub("viridis::scale_fill_viridis_c", "ggplot2::scale_fill_viridis_c", body_UO_fixed)
body(SFinx_UO) <- parse(text = body_UO_fixed)

body_BO <- deparse(body(SFinx_BO))
body_BO_fixed <- gsub("egg::background_image", "ggpubr::background_image", body_BO)
body_BO_fixed <- gsub("viridis::scale_fill_viridis_c", "ggplot2::scale_fill_viridis_c", body_BO_fixed)
body(SFinx_BO) <- parse(text = body_BO_fixed)

blank_image <- matrix("#FFFFFF", nrow = 1, ncol = 1)
```

### 2. Generate Synthetic Dataset

```r
n_spots <- 150

spatial_coords <- data.frame(
  Longitude = runif(n_spots, min = 0, max = 100),
  Latitude  = runif(n_spots, min = 0, max = 100)
)

normalized_scores <- data.frame(
  GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION = rnorm(n_spots, mean = 0, sd = 1.5),
  HP_LUPUS_NEPHRITIS                        = rnorm(n_spots, mean = 0, sd = 1.2)
)

combined_df <- cbind(spatial_coords, normalized_scores)
spatial_coords$x <- spatial_coords$Longitude
spatial_coords$y <- spatial_coords$Latitude
```

### 3. Test Unary Operator (SFinx_UO)

```r
u_plot <- SFinx_UO(
  Data   = normalized_scores,
  Finpre = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
  Thresh = 1.0,
  Pos    = spatial_coords,
  Image  = blank_image,
  Title  = "Neutrophil Activation Surface"
)
print(u_plot)
```

![Unary Operator Surface](https://github.com/mortezahaji/SFinxSDFA/blob/main/images/Unary_Operator_Plot.png)

### 4. Test Binary Operators (SFinx_BO)

```r
binary_res <- SFinx_BO(
  Data    = normalized_scores,
  Finpre1 = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
  Finpre2 = "HP_LUPUS_NEPHRITIS",
  Thresh  = 1.0,
  Pos     = spatial_coords,
  Image   = blank_image
)
print(binary_res$combined)
```

![Binary Operators Crosstalk](https://github.com/mortezahaji/SFinxSDFA/blob/main/images/Binary_Operators_Plot.png)

### 5. Fit Spatial Functional Regression (FOF_fit)

```r
grid_resolution <- 30

sdfa_fit <- FOF_fit(
  Data      = combined_df,
  Response  = "HP_LUPUS_NEPHRITIS",
  Predictor = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
  k         = 12,
  grid_n    = grid_resolution
)
```

### 6. Extract Surface Coefficients & Generate 3D Plots

```r
beta_matrix <- sdfa_fit$beta_surface
grid_x <- seq(min(combined_df$Longitude), max(combined_df$Longitude), length.out = grid_resolution)
grid_y <- seq(min(combined_df$Latitude),  max(combined_df$Latitude),  length.out = grid_resolution)

# Interactive 3D Plot
p_3d_interactive <- plot_ly(x = ~grid_x, y = ~grid_y, z = ~beta_matrix) %>%
  add_surface(colorscale = "Viridis",
              contours = list(z = list(show = TRUE, usecolormap = TRUE,
                                       highlightcolor = "#ff0000",
                                       project = list(z = TRUE)))) %>%
  layout(title = list(text = "3D Spatially Varying Coefficient Surface β(s)"),
         scene = list(xaxis = list(title = "Longitude"),
                      yaxis = list(title = "Latitude"),
                      zaxis = list(title = "Beta Coefficient Value")))
print(p_3d_interactive)
```

![FOF_fit 3D Surface](https://github.com/mortezahaji/SFinxSDFA/blob/main/images/3D_Coefficient_Surface.png)

---

## Citation

If you use SFinxSDFA in your research, please cite our preprint:

Hajihosseini, M., Patino-Martinez, E., Ghosal, R., Kaplan, M. J., & Pyne, S. (2026).
An automated platform for spatial functional modeling and fingerprint analysis of tissue
molecular landscapes. *bioRxiv* 2026.07.15.738836;
doi: [10.64898/2026.07.15.738836](https://doi.org/10.64898/2026.07.15.738836)
