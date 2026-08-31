# ==============================================================================
# Pure-R exact reproduction of Panels F and G
# Figure 3F-G analysis workflow
# This script ports the original Python calculations to R
# 1. Stable-feature filtering across all 12 treatment-year groups
# 2. Half-minimum imputation, sample-median normalization and log2 conversion
# 3. Original custom empirical-Bayes-like moderated F calculation
# 4. Year-adjusted Spearman associations
# 5. Original strict/supporting screening and 4 + 4 candidate selection
# 6. Original candidate-based temporal pathway-module analysis
# ==============================================================================
options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260717)
# ------------------------------------------------------------------------------
# 0. Package preparation / R
# ------------------------------------------------------------------------------
required_packages <- c("readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}
suppressPackageStartupMessages({
  library(readxl)
})
# ------------------------------------------------------------------------------
# 1. Input files / 4
# ------------------------------------------------------------------------------
soil_file <- "main_plant_soil_data.xlsx"
species_file <- "species_abundance.csv"
pos_file <- "metabolomics_positive.xlsx"
neg_file <- "metabolomics_negative.xlsx"
output_dir <- "pure_R_exact_FG_results-Final"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
required_files <- c(
  soil_file,
  species_file,
  pos_file,
  neg_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    " / Missing files: ",
    paste(missing_files, collapse = ", ")
  )
}
taxa <- c(
  "Staphylococcus capitis",
  "Pseudomonas syncyanea",
  "Luteitalea pratensis",
  "Paraflavitalea soli"
)
years <- c(0, 2, 4, 8, 12, 16)
positive_colour <- "#D65A7A"
negative_colour <- "#4F83D1"
high_header <- "#F5D9DE"
low_header <- "#DDEAF6"
canonical_sample <- function(x) {
  x <- sub("^RT-HCK-", "RT-H0-", x)
  x <- sub("^RT-GCK-", "RT-G0-", x)
  x
}
# ------------------------------------------------------------------------------
# 2. Read the four raw files / 4
# ------------------------------------------------------------------------------
plant_raw <- read_xlsx(
  soil_file,
  sheet = "Plant",
  .name_repair = "minimal"
)
species_raw <- read.csv(
  species_file,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
positive_raw <- read_xlsx(
  pos_file,
  sheet = "Sheet1",
  .name_repair = "minimal"
)
negative_raw <- read_xlsx(
  neg_file,
  sheet = "Sheet1",
  .name_repair = "minimal"
)
if (!"Species" %in% names(species_raw)) {
  stop(" Species.")
}
if (!all(taxa %in% species_raw$Species)) {
  stop(
    ":",
    paste(setdiff(taxa, species_raw$Species), collapse = ", ")
  )
}
# ------------------------------------------------------------------------------
# 3. Construct plant and bacterial metadata
# ------------------------------------------------------------------------------
plant_key <- as.character(plant_raw[["AnalysisID"]])
plant_lookup <- data.frame(
  Key = plant_key,
  Yield = as.numeric(plant_raw[["Yield_t_ha"]]),
  Treatment = substr(plant_key, 1, 1),
  Year = as.integer(plant_raw[["Year"]]),
  Rep = as.integer(sub(".*-", "", plant_key)),
  stringsAsFactors = FALSE
)
rownames(plant_lookup) <- plant_lookup$Key
species_sample_raw <- names(species_raw)[-1]
species_sample_canonical <- canonical_sample(species_sample_raw)
if (anyDuplicated(species_sample_canonical)) {
  stop("CK→0.")
}
species_column_by_sample <- setNames(
  species_sample_raw,
  species_sample_canonical
)
positive_sample_columns <- names(positive_raw)[25:ncol(positive_raw)]
negative_sample_columns <- names(negative_raw)[25:ncol(negative_raw)]
if (!identical(positive_sample_columns, negative_sample_columns)) {
  stop(",.")
}
canonical_metabolome_samples <- canonical_sample(
  positive_sample_columns
)
metadata_list <- vector(
  "list",
  length(canonical_metabolome_samples)
)
for (i in seq_along(canonical_metabolome_samples)) {
  sample_name <- canonical_metabolome_samples[i]
  key <- sub("^RT-", "", sample_name)
  if (!key %in% rownames(plant_lookup)) {
    stop(":", key)
  }
  if (!sample_name %in% names(species_column_by_sample)) {
    stop(":", sample_name)
  }
  species_column <- species_column_by_sample[[sample_name]]
  bacterial_values <- vapply(
    taxa,
    function(taxon) {
      row_index <- match(taxon, species_raw$Species)
      as.numeric(species_raw[row_index, species_column])
    },
    numeric(1)
  )
  plant_row <- plant_lookup[key, , drop = FALSE]
  metadata_list[[i]] <- data.frame(
    Sample = sample_name,
    Treatment = plant_row$Treatment,
    Year = plant_row$Year,
    Rep = plant_row$Rep,
    Yield = plant_row$Yield,
    Staphylococcus_capitis = bacterial_values[1],
    Pseudomonas_syncyanea = bacterial_values[2],
    Luteitalea_pratensis = bacterial_values[3],
    Paraflavitalea_soli = bacterial_values[4],
    stringsAsFactors = FALSE
  )
}
metadata <- do.call(rbind, metadata_list)
rownames(metadata) <- NULL
if (nrow(metadata) != 39) {
  stop(
    "39,:",
    nrow(metadata)
  )
}
# ------------------------------------------------------------------------------
# 4. Convert raw metabolomic tables to numeric matrices
# ------------------------------------------------------------------------------
to_numeric_matrix <- function(data) {
  output <- as.matrix(data)
  storage.mode(output) <- "numeric"
  output[!is.finite(output)] <- 0
  output
}
positive_annotations_all <- as.matrix(
  positive_raw[, 1:24, drop = FALSE]
)
negative_annotations_all <- as.matrix(
  negative_raw[, 1:24, drop = FALSE]
)
positive_matrix_all <- to_numeric_matrix(
  positive_raw[, 25:ncol(positive_raw), drop = FALSE]
)
negative_matrix_all <- to_numeric_matrix(
  negative_raw[, 25:ncol(negative_raw), drop = FALSE]
)
# ------------------------------------------------------------------------------
# 5. Stable-feature filtering / 12
# ------------------------------------------------------------------------------
group_indices <- list()
for (treatment in c("G", "H")) {
  for (year in years) {
    group_name <- paste0(treatment, year)
    indices <- which(
      metadata$Treatment == treatment &
      metadata$Year == year
    )
    if (length(indices) == 0) {
      stop(":", group_name)
    }
    group_indices[[group_name]] <- indices
  }
}
stable_feature_filter <- function(matrix, groups) {
  detected <- matrix > 0
  keep <- rep(TRUE, nrow(matrix))
  for (indices in groups) {
    required <- ceiling(length(indices) * 0.50)
    keep <- keep & (
      rowSums(
        detected[, indices, drop = FALSE]
      ) >= required
    )
  }
  keep
}
positive_keep <- stable_feature_filter(
  positive_matrix_all,
  group_indices
)
negative_keep <- stable_feature_filter(
  negative_matrix_all,
  group_indices
)
if (sum(positive_keep) != 67) {
  stop(
    "POS67,:",
    sum(positive_keep)
  )
}
if (sum(negative_keep) != 28) {
  stop(
    "NEG28,:",
    sum(negative_keep)
  )
}
# ------------------------------------------------------------------------------
# 6. Imputation, median normalization and log2 transformation
# ,log2
# ------------------------------------------------------------------------------
impute_normalize_log2 <- function(matrix) {
  output <- matrix
  for (row_index in seq_len(nrow(output))) {
    nonzero <- output[row_index, output[row_index, ] > 0]
    replacement <- if (length(nonzero) > 0) {
      min(nonzero) / 2
    } else {
      1
    }
    output[
      row_index,
      output[row_index, ] <= 0
    ] <- replacement
  }
  sample_medians <- apply(output, 2, median)
  reference_median <- median(sample_medians)
  normalized <- sweep(
    output,
    2,
    sample_medians,
    "/"
  ) * reference_median
  log2(normalized)
}
positive_annotations <- positive_annotations_all[
  positive_keep,
  ,
  drop = FALSE
]
negative_annotations <- negative_annotations_all[
  negative_keep,
  ,
  drop = FALSE
]
positive_log <- impute_normalize_log2(
  positive_matrix_all[
    positive_keep,
    ,
    drop = FALSE
  ]
)
negative_log <- impute_normalize_log2(
  negative_matrix_all[
    negative_keep,
    ,
    drop = FALSE
  ]
)
# samples × features
metabolite_matrix <- t(
  rbind(
    positive_log,
    negative_log
  )
)
annotations <- rbind(
  positive_annotations,
  negative_annotations
)
modes <- c(
  rep("POS", nrow(positive_annotations)),
  rep("NEG", nrow(negative_annotations))
)
# ------------------------------------------------------------------------------
# 7. Original custom empirical-Bayes-like mowing model
# ------------------------------------------------------------------------------
group_order <- list(
  c("G", 0), c("H", 0),
  c("G", 2), c("H", 2),
  c("G", 4), c("H", 4),
  c("G", 8), c("H", 8),
  c("G", 12), c("H", 12),
  c("G", 16), c("H", 16)
)
design <- sapply(
  group_order,
  function(group) {
    as.numeric(
      metadata$Treatment == group[1] &
      metadata$Year == as.numeric(group[2])
    )
  }
)
xtx_inverse <- solve(crossprod(design))
coefficients <- (
  xtx_inverse %*%
  t(design) %*%
  metabolite_matrix
)
residuals_matrix <- (
  metabolite_matrix -
  design %*% coefficients
)
residual_df <- (
  nrow(metabolite_matrix) -
  qr(design)$rank
)
raw_variance <- (
  colSums(residuals_matrix^2) /
  residual_df
)
post_contrasts <- matrix(
  0,
  nrow = 12,
  ncol = 5
)
contrast_columns <- list(
  c(3, 4),
  c(5, 6),
  c(7, 8),
  c(9, 10),
  c(11, 12)
)
for (contrast_index in seq_along(contrast_columns)) {
  g_column <- contrast_columns[[contrast_index]][1]
  h_column <- contrast_columns[[contrast_index]][2]
  post_contrasts[
    h_column,
    contrast_index
  ] <- 1
  post_contrasts[
    g_column,
    contrast_index
  ] <- -1
}
post_effects <- (
  t(post_contrasts) %*%
  coefficients
)
contrast_covariance <- (
  t(post_contrasts) %*%
  xtx_inverse %*%
  post_contrasts
)
contrast_covariance_inverse <- solve(
  contrast_covariance
)
quadratic_form <- colSums(
  post_effects *
  (
    contrast_covariance_inverse %*%
    post_effects
  )
)
squeeze_variances <- function(
  variances,
  df
) {
  adjusted_log_variance <- (
    log(variances) -
    digamma(df / 2) +
    log(df / 2)
  )
  target <- (
    var(adjusted_log_variance) -
    trigamma(df / 2)
  )
  if (target <= 1e-8) {
    prior_df <- 1e6
  } else {
    prior_df <- uniroot(
      function(value) {
        trigamma(value / 2) - target
      },
      interval = c(0.01, 1e6)
    )$root
  }
  prior_variance <- exp(
    mean(adjusted_log_variance) +
    digamma(prior_df / 2) -
    log(prior_df / 2)
  )
  posterior <- (
    df * variances +
    prior_df * prior_variance
  ) / (
    df + prior_df
  )
  list(
    posterior = posterior,
    prior_df = prior_df
  )
}
variance_fit <- squeeze_variances(
  raw_variance,
  residual_df
)
posterior_variance <- variance_fit$posterior
prior_df <- variance_fit$prior_df
moderated_f <- (
  quadratic_form / 5
) / posterior_variance
mowing_p <- pf(
  moderated_f,
  df1 = 5,
  df2 = residual_df + prior_df,
  lower.tail = FALSE
)
mowing_fdr <- p.adjust(
  mowing_p,
  method = "BH"
)
mean_mowing_effect <- colMeans(post_effects)
dominant_sign <- sign(
  apply(post_effects, 2, median)
)
dominant_sign[dominant_sign == 0] <- 1
same_direction_years <- colSums(
  sweep(
    post_effects,
    2,
    dominant_sign,
    "*"
  ) > 0
)
fold_years <- colSums(
  abs(post_effects) >= log2(1.2)
)
strong_mowing <- (
  mowing_fdr < 0.05 &
  same_direction_years >= 4 &
  fold_years >= 2
)
# ------------------------------------------------------------------------------
# 8. Original year-adjusted Spearman associations
# Spearman
# ------------------------------------------------------------------------------
year_design <- cbind(
  Intercept = 1,
  sapply(
    years[-1],
    function(current_year) {
      as.numeric(metadata$Year == current_year)
    }
  )
)
year_projection <- (
  year_design %*%
  solve(crossprod(year_design)) %*%
  t(year_design)
)
year_adjusted_spearman_matrix <- function(
  matrix,
  target
) {
  rank_matrix <- apply(
    matrix,
    2,
    rank,
    ties.method = "average"
  )
  if (is.null(dim(rank_matrix))) {
    rank_matrix <- matrix(
      rank_matrix,
      ncol = 1
    )
  }
  rank_target <- rank(
    target,
    ties.method = "average"
  )
  matrix_residual <- (
    rank_matrix -
    year_projection %*% rank_matrix
  )
  target_residual <- (
    rank_target -
    year_projection %*% rank_target
  )
  numerator <- colSums(
    sweep(
      matrix_residual,
      1,
      target_residual,
      "*"
    )
  )
  denominator <- sqrt(
    colSums(matrix_residual^2) *
    sum(target_residual^2)
  )
  rho <- ifelse(
    denominator > 0,
    numerator / denominator,
    0
  )
  df <- (
    nrow(matrix) -
    qr(year_design)$rank -
    1
  )
  t_value <- rho * sqrt(
    df /
    pmax(1 - rho^2, 1e-12)
  )
  p_value <- 2 * pt(
    abs(t_value),
    df = df,
    lower.tail = FALSE
  )
  q_value <- p.adjust(
    p_value,
    method = "BH"
  )
  list(
    rho = rho,
    p = p_value,
    fdr = q_value
  )
}
yield_vector <- metadata$Yield
bacteria_matrix <- log1p(
  as.matrix(
    metadata[, c(
      "Staphylococcus_capitis",
      "Pseudomonas_syncyanea",
      "Luteitalea_pratensis",
      "Paraflavitalea_soli"
    )]
  )
)
bacteria_z <- scale(
  bacteria_matrix,
  center = TRUE,
  scale = TRUE
)
bacterial_axis <- (
  bacteria_z[, 1] +
  bacteria_z[, 2] -
  bacteria_z[, 3] -
  bacteria_z[, 4]
) / 2
yield_association <- year_adjusted_spearman_matrix(
  metabolite_matrix,
  yield_vector
)
axis_association <- year_adjusted_spearman_matrix(
  metabolite_matrix,
  bacterial_axis
)
taxon_associations <- lapply(
  seq_len(4),
  function(taxon_index) {
    year_adjusted_spearman_matrix(
      metabolite_matrix,
      bacteria_matrix[, taxon_index]
    )
  }
)
yield_rho <- yield_association$rho
yield_fdr <- yield_association$fdr
axis_rho <- axis_association$rho
axis_fdr <- axis_association$fdr
taxon_rho <- do.call(
  rbind,
  lapply(taxon_associations, `[[`, "rho")
)
taxon_fdr <- do.call(
  rbind,
  lapply(taxon_associations, `[[`, "fdr")
)
direction_aligned <- (
  sign(mean_mowing_effect) ==
  sign(yield_rho)
) & (
  sign(mean_mowing_effect) ==
  sign(axis_rho)
)
expected_signs <- rbind(
  sign(yield_rho),
  sign(yield_rho),
  -sign(yield_rho),
  -sign(yield_rho)
)
bacteria_direction_agreement <- colSums(
  sign(taxon_rho) == expected_signs
)
strict_core <- (
  strong_mowing &
  yield_fdr < 0.05 &
  axis_fdr < 0.05 &
  direction_aligned &
  bacteria_direction_agreement >= 3
)
supporting <- (
  strong_mowing &
  yield_fdr < 0.10 &
  axis_fdr < 0.10 &
  direction_aligned &
  bacteria_direction_agreement >= 3 &
  !strict_core
)
# ------------------------------------------------------------------------------
# 9. Build the metabolite-result table
# ------------------------------------------------------------------------------
annotation_fields <- function(
  annotation,
  mode
) {
  full_match_columns <- annotation[18:20]
  data.frame(
    FeatureID = as.character(annotation[1]),
    Name = as.character(annotation[2]),
    Mode = mode,
    KEGG = as.character(annotation[8]),
    HMDB = as.character(annotation[10]),
    Annotation_full_matches = sum(
      grepl(
        "Full match",
        as.character(full_match_columns),
        fixed = TRUE
      ),
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}
metabolite_results_list <- vector(
  "list",
  ncol(metabolite_matrix)
)
for (feature_index in seq_len(ncol(metabolite_matrix))) {
  base <- annotation_fields(
    annotations[feature_index, ],
    modes[feature_index]
  )
  mowing_fdr_value <- mowing_fdr[feature_index]
  joint_score <- (
    min(
      -log10(
        max(mowing_fdr_value, 1e-300)
      ),
      10
    ) *
    (1 + abs(mean_mowing_effect[feature_index])) *
    (1 + abs(yield_rho[feature_index])) *
    (1 + abs(axis_rho[feature_index])) *
    (same_direction_years[feature_index] / 5)
  )
  metabolite_results_list[[feature_index]] <- data.frame(
    base,
    FeatureIndex = feature_index,
    Mowing_effect = mean_mowing_effect[feature_index],
    Mowing_FDR = mowing_fdr_value,
    Yield_rho = yield_rho[feature_index],
    Yield_FDR = yield_fdr[feature_index],
    Axis_rho = axis_rho[feature_index],
    Axis_FDR = axis_fdr[feature_index],
    S_capitis_rho = taxon_rho[1, feature_index],
    P_syncyanea_rho = taxon_rho[2, feature_index],
    L_pratensis_rho = taxon_rho[3, feature_index],
    P_soli_rho = taxon_rho[4, feature_index],
    Strict_core = strict_core[feature_index],
    Supporting = supporting[feature_index],
    Direction_agreement_n = bacteria_direction_agreement[feature_index],
    Same_direction_years = same_direction_years[feature_index],
    Strong_mowing = strong_mowing[feature_index],
    Joint_score = joint_score,
    stringsAsFactors = FALSE
  )
}
metabolite_results <- do.call(
  rbind,
  metabolite_results_list
)
# ------------------------------------------------------------------------------
# 10. Original candidate selection, including 4 + 4 completion
# Select four high-yield and four low-yield candidates
# ------------------------------------------------------------------------------
candidate_pool <- metabolite_results[
  metabolite_results$Strict_core |
  metabolite_results$Supporting,
  ,
  drop = FALSE
]
candidate_pool <- candidate_pool[
  order(-candidate_pool$Joint_score),
  ,
  drop = FALSE
]
candidate_pool <- candidate_pool[
  !duplicated(candidate_pool$Name),
  ,
  drop = FALSE
]
high_pool <- candidate_pool[
  candidate_pool$Mowing_effect > 0 &
  candidate_pool$Yield_rho > 0 &
  candidate_pool$Axis_rho > 0,
  ,
  drop = FALSE
]
high_pool <- high_pool[
  order(
    !high_pool$Strict_core,
    -high_pool$Joint_score
  ),
  ,
  drop = FALSE
]
high_candidates <- head(high_pool, 4)
low_pool <- candidate_pool[
  candidate_pool$Mowing_effect < 0 &
  candidate_pool$Yield_rho < 0 &
  candidate_pool$Axis_rho < 0,
  ,
  drop = FALSE
]
low_pool <- low_pool[
  order(
    !low_pool$Strict_core,
    -low_pool$Joint_score
  ),
  ,
  drop = FALSE
]
low_candidates <- head(low_pool, 4)
fill_candidates <- function(
  selected,
  positive
) {
  selected_names <- selected$Name
  if (positive) {
    alternatives <- metabolite_results[
      metabolite_results$Strong_mowing &
      metabolite_results$Mowing_effect > 0 &
      metabolite_results$Yield_rho > 0 &
      metabolite_results$Axis_rho > 0 &
      !metabolite_results$Name %in% selected_names,
      ,
      drop = FALSE
    ]
  } else {
    alternatives <- metabolite_results[
      metabolite_results$Strong_mowing &
      metabolite_results$Mowing_effect < 0 &
      metabolite_results$Yield_rho < 0 &
      metabolite_results$Axis_rho < 0 &
      !metabolite_results$Name %in% selected_names,
      ,
      drop = FALSE
    ]
  }
  alternatives <- alternatives[
    order(-alternatives$Joint_score),
    ,
    drop = FALSE
  ]
  alternatives <- alternatives[
    !duplicated(alternatives$Name),
    ,
    drop = FALSE
  ]
  while (
    nrow(selected) < 4 &&
    nrow(alternatives) > 0
  ) {
    selected <- rbind(
      selected,
      alternatives[1, , drop = FALSE]
    )
    alternatives <- alternatives[
      alternatives$Name != selected$Name[nrow(selected)],
      ,
      drop = FALSE
    ]
  }
  head(selected, 4)
}
high_candidates <- fill_candidates(
  high_candidates,
  TRUE
)
low_candidates <- fill_candidates(
  low_candidates,
  FALSE
)
figure_f_rows <- rbind(
  high_candidates,
  low_candidates
)
# ------------------------------------------------------------------------------
# 11. Mandatory verification
# ------------------------------------------------------------------------------
expected_high <- c(
  "Sakuranetin",
  "Hypoxanthine",
  "Asp-Glu",
  "L-Kynurenine"
)
expected_low <- c(
  "Uridine",
  "Catalpol",
  "Tyrosol",
  "8-O-Acetylharpagide"
)
actual_high <- high_candidates$Name
actual_low <- low_candidates$Name
verification <- data.frame(
  Item = c(
    "Metabolomic samples",
    "Stable POS features",
    "Stable NEG features",
    "Strong mowing features",
    "Strict core features",
    "Supporting features"
  ),
  Value = c(
    nrow(metadata),
    sum(positive_keep),
    sum(negative_keep),
    sum(strong_mowing),
    sum(strict_core),
    sum(supporting)
  ),
  Expected = c(
    39,
    67,
    28,
    29,
    6,
    3
  )
)
if (!all(verification$Value == verification$Expected)) {
  print(verification)
  stop(",.")
}
if (!identical(actual_high, expected_high)) {
  stop(
    ".\nExpected: ",
    paste(expected_high, collapse = ", "),
    "\nActual: ",
    paste(actual_high, collapse = ", ")
  )
}
if (!identical(actual_low, expected_low)) {
  stop(
    ".\nExpected: ",
    paste(expected_low, collapse = ", "),
    "\nActual: ",
    paste(actual_low, collapse = ", ")
  )
}
write.csv(
  metabolite_results,
  file.path(
    output_dir,
    "All_metabolite_screening_results.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  figure_f_rows,
  file.path(
    output_dir,
    "Figure_F_candidate_metabolites.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  verification,
  file.path(
    output_dir,
    "verification.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
# ------------------------------------------------------------------------------
# 12. Panel F drawing function / F
# ------------------------------------------------------------------------------
draw_panel_f <- function(
  high_candidates,
  low_candidates,
  panel_letter = "F"
) {
  value_columns <- c(
    "Mowing_effect",
    "Yield_rho",
    "S_capitis_rho",
    "P_syncyanea_rho",
    "L_pratensis_rho",
    "P_soli_rho"
  )
  column_labels <- list(
    "Mowing",
    "Yield",
    expression(italic(S.)~italic(capitis)),
    expression(italic(P.)~italic(syncyanea)),
    expression(italic(L.)~italic(pratensis)),
    expression(italic(P.)~italic(soli))
  )
  n_high <- nrow(high_candidates)
  n_low <- nrow(low_candidates)
  left_edge <- 0.10
  label_right <- 2.75
  matrix_left <- 2.75
  cell_width <- 0.78
  matrix_right <- matrix_left + cell_width * length(value_columns)
  title_y <- 10.92
  column_y <- 10.18
  table_top <- 9.82
  header_height <- 0.72
  row_height <- 0.72
  high_header_top <- table_top
  high_header_bottom <- high_header_top - header_height
  high_rows_top <- high_header_bottom
  high_rows_bottom <- high_rows_top - n_high * row_height
  low_header_top <- high_rows_bottom
  low_header_bottom <- low_header_top - header_height
  low_rows_top <- low_header_bottom
  low_rows_bottom <- low_rows_top - n_low * row_height
  legend_y <- low_rows_bottom - 0.62
  plot.new()
  plot.window(
    xlim = c(-0.25, matrix_right + 0.12),
    ylim = c(legend_y - 0.35, title_y + 0.28),
    xaxs = "i",
    yaxs = "i"
  )
  text(
    -0.18,
    title_y,
    panel_letter,
    cex = 1.14,
    font = 2,
    adj = c(0, 0.5)
  )
  text(
    (left_edge + matrix_right) / 2,
    title_y,
    "Candidate root-tissue metabolites",
    cex = 0.92,
    font = 2
  )
  for (j in seq_along(column_labels)) {
    x <- matrix_left + cell_width * (j - 0.5)
    text(
      x,
      column_y,
      labels = column_labels[[j]],
      cex = 0.67,
      adj = c(0.5, 0)
    )
  }
  rect(
    left_edge,
    high_header_bottom,
    matrix_right,
    high_header_top,
    col = high_header,
    border = NA
  )
  text(
    (left_edge + matrix_right) / 2,
    (high_header_top + high_header_bottom) / 2,
    "High-yield associated (mowing-enriched)",
    cex = 0.70
  )
  rect(
    left_edge,
    low_header_bottom,
    matrix_right,
    low_header_top,
    col = low_header,
    border = NA
  )
  text(
    (left_edge + matrix_right) / 2,
    (low_header_top + low_header_bottom) / 2,
    "Low-yield associated (mowing-reduced)",
    cex = 0.70
  )
  draw_rows <- function(data, top_y) {
    for (i in seq_len(nrow(data))) {
      y <- top_y - row_height * (i - 0.5)
      display_name <- as.character(data$Name[i])
      if (
        "Supporting" %in% names(data) &&
        isTRUE(data$Supporting[i])
      ) {
        display_name <- paste0(display_name, "†")
      }
      text(
        left_edge + 0.10,
        y,
        paste0(i, ". ", display_name),
        adj = c(0, 0.5),
        cex = 0.68,
        font = 1
      )
      values <- as.numeric(data[i, value_columns, drop = TRUE])
      for (j in seq_along(values)) {
        x <- matrix_left + cell_width * (j - 0.5)
        fill_col <- if (values[j] > 0) positive_colour else negative_colour
        points(
          x,
          y,
          pch = 21,
          bg = fill_col,
          col = "white",
          lwd = 0.55,
          cex = 1.15
        )
      }
    }
  }
  draw_rows(high_candidates, high_rows_top)
  draw_rows(low_candidates, low_rows_top)
  horizontal_lines <- c(
    table_top,
    high_header_bottom,
    high_rows_top - row_height * (0:n_high),
    low_header_bottom,
    low_rows_top - row_height * (0:n_low)
  )
  for (y in sort(unique(round(horizontal_lines, 8)), decreasing = TRUE)) {
    segments(
      left_edge,
      y,
      matrix_right,
      y,
      col = "grey72",
      lwd = 0.65
    )
  }
  vertical_lines <- c(
    left_edge,
    label_right,
    matrix_left + cell_width * (0:length(value_columns))
  )
  for (x in unique(round(vertical_lines, 8))) {
    segments(
      x,
      low_rows_bottom,
      x,
      table_top,
      col = "grey72",
      lwd = 0.65
    )
  }
  rect(
    left_edge,
    low_rows_bottom,
    matrix_right,
    table_top,
    border = "grey45",
    lwd = 0.8
  )
  points(
    1.18,
    legend_y,
    pch = 21,
    bg = positive_colour,
    col = "white",
    lwd = 0.55,
    cex = 1.00
  )
  text(
    1.40,
    legend_y,
    "Positive / increased",
    adj = c(0, 0.5),
    cex = 0.64
  )
  points(
    4.05,
    legend_y,
    pch = 21,
    bg = negative_colour,
    col = "white",
    lwd = 0.55,
    cex = 1.00
  )
  text(
    4.27,
    legend_y,
    "Negative / decreased",
    adj = c(0, 0.5),
    cex = 0.64
  )
  text(
    matrix_right,
    legend_y,
    "† Supporting candidate",
    adj = c(1, 0.5),
    cex = 0.61
  )
}
# ------------------------------------------------------------------------------
# 13. Original candidate-based yearwise pathway analysis for Panel G
# G
# ------------------------------------------------------------------------------
module_map <- list(
  "Flavonoid and phenolic metabolism" = c(
    "Sakuranetin",
    "Tyrosol"
  ),
  "Purine metabolism" = c(
    "Hypoxanthine"
  ),
  "Peptide and amino-acid metabolism" = c(
    "Asp-Glu"
  ),
  "Tryptophan–kynurenine metabolism" = c(
    "L-Kynurenine"
  ),
  "Pyrimidine metabolism" = c(
    "Uridine"
  ),
  "Iridoid glycoside metabolism" = c(
    "Catalpol",
    "8-O-Acetylharpagide"
  )
)
candidate_feature_values <- vector(
  "list",
  nrow(figure_f_rows)
)
for (row_index in seq_len(nrow(figure_f_rows))) {
  feature_id <- as.character(
    figure_f_rows$FeatureID[row_index]
  )
  mode <- figure_f_rows$Mode[row_index]
  if (mode == "POS") {
    source_data <- positive_raw
    id_column <- "ID"
  } else {
    source_data <- negative_raw
    id_column <- "Compound_ID"
  }
  matched_indices <- which(
    as.character(source_data[[id_column]]) ==
    feature_id
  )
  if (length(matched_indices) != 1) {
    stop(
      "FeatureID:",
      feature_id
    )
  }
  values <- as.numeric(
    source_data[
      matched_indices,
      positive_sample_columns,
      drop = TRUE
    ]
  )
  values[!is.finite(values)] <- 0
  nonzero <- values[values > 0]
  replacement <- if (
    length(nonzero) > 0
  ) {
    min(nonzero) / 2
  } else {
    1
  }
  values[values <= 0] <- replacement
  candidate_feature_values[[row_index]] <- log2(
    values
  )
}
yearwise_list <- list()
counter <- 1
for (row_index in seq_len(nrow(figure_f_rows))) {
  metabolite_name <- figure_f_rows$Name[row_index]
  feature_id <- figure_f_rows$FeatureID[row_index]
  values <- candidate_feature_values[[row_index]]
  pathway_name <- names(
    Filter(
      function(metabolites) {
        metabolite_name %in% metabolites
      },
      module_map
    )
  )
  if (length(pathway_name) != 1) {
    stop(
      ":",
      metabolite_name
    )
  }
  for (current_year in years) {
    h_indices <- which(
      metadata$Treatment == "H" &
      metadata$Year == current_year
    )
    g_indices <- which(
      metadata$Treatment == "G" &
      metadata$Year == current_year
    )
    h_values <- values[h_indices]
    g_values <- values[g_indices]
    effect <- (
      mean(h_values) -
      mean(g_values)
    )
    p_value <- tryCatch(
      t.test(
        h_values,
        g_values,
        var.equal = FALSE
      )$p.value,
      error = function(e) 1
    )
    yearwise_list[[counter]] <- data.frame(
      FeatureID = feature_id,
      Name = metabolite_name,
      Pathway = pathway_name,
      Year = current_year,
      log2FC_H_vs_G = effect,
      P = p_value,
      stringsAsFactors = FALSE
    )
    counter <- counter + 1
  }
}
yearwise_results <- do.call(
  rbind,
  yearwise_list
)
yearwise_results$FDR <- NA_real_
for (current_year in years) {
  indices <- which(
    yearwise_results$Year ==
    current_year
  )
  yearwise_results$FDR[indices] <- p.adjust(
    yearwise_results$P[indices],
    method = "BH"
  )
}
module_year_list <- list()
counter <- 1
for (pathway_name in names(module_map)) {
  for (current_year in years) {
    subset_data <- yearwise_results[
      yearwise_results$Pathway ==
        pathway_name &
      yearwise_results$Year ==
        current_year,
      ,
      drop = FALSE
    ]
    fisher_statistic <- -2 * sum(
      log(
        pmax(
          subset_data$P,
          1e-300
        )
      )
    )
    fisher_p <- pchisq(
      fisher_statistic,
      df = 2 * nrow(subset_data),
      lower.tail = FALSE
    )
    responsive_n <- sum(
      abs(
        subset_data$log2FC_H_vs_G
      ) >= log2(1.2) &
      subset_data$FDR < 0.10
    )
    module_year_list[[counter]] <- data.frame(
      Pathway = pathway_name,
      Year = current_year,
      Mean_log2FC = mean(
        subset_data$log2FC_H_vs_G
      ),
      Responsive_n = responsive_n,
      Combined_P = fisher_p,
      Metabolites = paste(
        module_map[[pathway_name]],
        collapse = ", "
      ),
      stringsAsFactors = FALSE
    )
    counter <- counter + 1
  }
}
module_year_results <- do.call(
  rbind,
  module_year_list
)
module_year_results$Combined_FDR <- p.adjust(
  module_year_results$Combined_P,
  method = "BH"
)
write.csv(
  yearwise_results,
  file.path(
    output_dir,
    "Figure_G_metabolite_year_details.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  module_year_results,
  file.path(
    output_dir,
    "Figure_G_module_year_summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
# ------------------------------------------------------------------------------
# 14. Panel G drawing function / G
# ------------------------------------------------------------------------------
draw_panel_g <- function(
  module_year_results,
  module_map,
  panel_letter = "G"
) {
  pathway_order <- c(
    "Flavonoid and phenolic metabolism",
    "Purine metabolism",
    "Peptide and amino-acid metabolism",
    "Tryptophan–kynurenine metabolism",
    "Pyrimidine metabolism",
    "Iridoid glycoside metabolism"
  )
  pathway_labels <- c(
    "Flavonoid and\nphenolic metabolism",
    "Purine\nmetabolism",
    "Peptide and amino-acid\nmetabolism",
    "Tryptophan-kynurenine\nmetabolism",
    "Pyrimidine\nmetabolism",
    "Iridoid glycoside\nmetabolism"
  )
  label_by_pathway <- setNames(
    pathway_labels,
    pathway_order
  )
  y_order <- rev(pathway_order)
  y_positions <- setNames(
    seq_along(y_order) - 1,
    y_order
  )
  x_positions <- setNames(
    seq_along(years) - 1,
    years
  )
  format_metabolite_label <- function(x) {
    if (grepl(", ", x, fixed = TRUE) && nchar(x) > 18) {
      return(gsub(", ", ",\n", x, fixed = TRUE))
    }
    x
  }
  plot.new()
  plot.window(
    xlim = c(-0.80, 12.45),
    ylim = c(-1.95, length(y_order) - 0.40),
    xaxs = "i",
    yaxs = "i"
  )
  for (x in seq_along(years) - 1) {
    abline(
      v = x,
      col = "grey92",
      lwd = 0.7
    )
  }
  for (y in seq_along(y_order) - 1) {
    abline(
      h = y,
      col = "grey92",
      lwd = 0.7
    )
  }
  value_range <- range(
    module_year_results$Mean_log2FC,
    finite = TRUE
  )
  palette <- colorRampPalette(
    c(
      "#3B4CC0",
      "#FFFFFF",
      "#B40426"
    )
  )(256)
  colour_for_value <- function(value) {
    if (
      !is.finite(value) ||
      diff(value_range) == 0
    ) {
      return(palette[128])
    }
    scaled <- (
      value - value_range[1]
    ) / diff(value_range)
    index <- round(
      1 + scaled * 255
    )
    index <- min(
      max(index, 1),
      256
    )
    palette[index]
  }
  for (
    row_index in seq_len(
      nrow(module_year_results)
    )
  ) {
    row <- module_year_results[
      row_index,
      ,
      drop = FALSE
    ]
    x <- x_positions[
      as.character(row$Year)
    ]
    y <- y_positions[
      row$Pathway
    ]
    point_cex <- (
      1.10 +
      0.90 * row$Responsive_n
    )
    points(
      x,
      y,
      pch = 21,
      bg = colour_for_value(
        row$Mean_log2FC
      ),
      col = "black",
      lwd = 0.7,
      cex = point_cex
    )
    significance <- if (
      row$Combined_FDR < 0.01
    ) {
      "**"
    } else if (
      row$Combined_FDR < 0.05
    ) {
      "*"
    } else {
      ""
    }
    if (nzchar(significance)) {
      text(
        x,
        y,
        significance,
        cex = 0.74,
        font = 2
      )
    }
  }
  axis(
    1,
    at = seq_along(years) - 1,
    labels = years,
    cex.axis = 0.72,
    lwd = 0.7,
    tck = -0.02
  )
  axis(
    2,
    at = seq_along(y_order) - 1,
    labels = label_by_pathway[y_order],
    las = 1,
    cex.axis = 0.60,
    lwd = 0.7,
    tck = -0.02
  )
  mtext(
    "Cultivation year",
    side = 1,
    line = 2.2,
    cex = 0.76
  )
  mtext(
    "Pathway module",
    side = 2,
    line = 6.1,
    cex = 0.76
  )
  text(
    7.75,
    length(y_order) - 0.28,
    "Representative metabolites",
    adj = c(0, 0),
    cex = 0.67,
    font = 2
  )
  for (pathway_name in pathway_order) {
    met_label <- format_metabolite_label(
      paste(
        module_map[[pathway_name]],
        collapse = ", "
      )
    )
    text(
      7.75,
      y_positions[pathway_name],
      met_label,
      adj = c(0, 0.5),
      cex = 0.58,
      font = 3
    )
  }
 # Bubble-size legend placed below the matrix and kept clear of the y labels
  legend_counts <- c(0, 1, 2)
  legend_x <- c(1.05, 2.75, 4.45)
  legend_y <- -1.32
  text(
    -0.58,
    legend_y,
    "Responsive metabolites",
    adj = c(0, 0.5),
    cex = 0.56,
    font = 2
  )
  for (i in seq_along(legend_counts)) {
    count <- legend_counts[i]
    points(
      legend_x[i],
      legend_y,
      pch = 21,
      bg = "grey82",
      col = "black",
      cex = 1.10 + 0.90 * count
    )
    text(
      legend_x[i] + 0.28,
      legend_y,
      as.character(count),
      adj = c(0, 0.5),
      cex = 0.56
    )
  }
 # Colour bar moved further right to avoid overlap
  bar_x1 <- 10.85
  bar_x2 <- 11.10
  bar_y1 <- 0.45
  bar_y2 <- 4.55
  bar_breaks <- seq(
    bar_y1,
    bar_y2,
    length.out = 257
  )
  for (i in seq_len(256)) {
    rect(
      bar_x1,
      bar_breaks[i],
      bar_x2,
      bar_breaks[i + 1],
      col = palette[i],
      border = NA
    )
  }
  rect(
    bar_x1,
    bar_y1,
    bar_x2,
    bar_y2,
    border = "black",
    lwd = 0.6
  )
  tick_values <- pretty(
    value_range,
    n = 5
  )
  tick_y <- (
    bar_y1 +
    (
      tick_values - value_range[1]
    ) / diff(value_range) *
    (bar_y2 - bar_y1)
  )
  segments(
    bar_x2,
    tick_y,
    bar_x2 + 0.10,
    tick_y,
    lwd = 0.6
  )
  text(
    bar_x2 + 0.15,
    tick_y,
    format(
      tick_values,
      digits = 2,
      trim = TRUE
    ),
    adj = c(0, 0.5),
    cex = 0.54
  )
  text(
    11.80,
    (bar_y1 + bar_y2) / 2,
    "Mean mowing effect,\nlog2(H/G)",
    srt = 90,
    cex = 0.57
  )
  text(
    -0.66,
    length(y_order) - 0.06,
    panel_letter,
    cex = 1.15,
    font = 2,
    adj = c(0, 0.5)
  )
  text(
    12.28,
    -1.72,
    "* FDR < 0.05; ** FDR < 0.01",
    adj = c(1, 0.5),
    cex = 0.50
  )
}
# ------------------------------------------------------------------------------
# 15. Export Panel F, Panel G and combined F-G figure
# F,G
# ------------------------------------------------------------------------------
export_plot <- function(
  filename,
  width,
  height,
  plot_function,
  type,
  mar = c(0, 0, 0, 0)
) {
  if (type == "svg") {
    svg(
      filename,
      width = width,
      height = height,
      family = "sans",
      bg = "transparent"
    )
  } else if (type == "png") {
    png(
      filename,
      width = width,
      height = height,
      units = "in",
      res = 600,
      bg = "transparent"
    )
  } else if (type == "tiff") {
    tiff(
      filename,
      width = width,
      height = height,
      units = "in",
      res = 600,
      compression = "lzw",
      bg = "transparent"
    )
  } else {
    stop("Unknown output type.")
  }
  par(
    mar = mar,
    family = "sans",
    xpd = NA
  )
  plot_function()
  dev.off()
}
for (type in c("svg", "png", "tiff")) {
  export_plot(
    file.path(
      output_dir,
      paste0(
        "FigF_candidate_metabolite_matrix.",
        type
      )
    ),
    width = 7.20,
    height = 4.90,
    plot_function = function() {
      draw_panel_f(
        high_candidates,
        low_candidates,
        "F"
      )
    },
    type = type,
    mar = c(0, 0, 0, 0)
  )
  export_plot(
    file.path(
      output_dir,
      paste0(
        "FigG_pathway_module_temporal_response.",
        type
      )
    ),
    width = 8.90,
    height = 6.20,
    plot_function = function() {
      draw_panel_g(
        module_year_results,
        module_map,
        "G"
      )
    },
    type = type,
    mar = c(2.9, 7.8, 1.1, 2.1)
  )
}
# Combined F-G
for (type in c("svg", "png", "tiff")) {
  filename <- file.path(
    output_dir,
    paste0(
      "FigFG_combined.",
      type
    )
  )
  if (type == "svg") {
    svg(
      filename,
      width = 18.20,
      height = 6.00,
      family = "sans",
      bg = "transparent"
    )
  } else if (type == "png") {
    png(
      filename,
      width = 18.20,
      height = 6.00,
      units = "in",
      res = 600,
      bg = "transparent"
    )
  } else {
    tiff(
      filename,
      width = 18.20,
      height = 6.00,
      units = "in",
      res = 600,
      compression = "lzw",
      bg = "transparent"
    )
  }
  layout(
    matrix(c(1, 2, 3), nrow = 1),
    widths = c(0.95, 0.08, 1.48)
  )
  par(
    mar = c(0, 0, 0, 0),
    family = "sans",
    xpd = NA
  )
  draw_panel_f(
    high_candidates,
    low_candidates,
    "F"
  )
 # spacer panel to prevent overlap between F and G
 # ,FG
  par(
    mar = c(0, 0, 0, 0),
    family = "sans",
    xpd = NA
  )
  plot.new()
  par(
    mar = c(2.7, 7.6, 1.0, 2.0),
    family = "sans",
    xpd = NA
  )
  draw_panel_g(
    module_year_results,
    module_map,
    "G"
  )
  dev.off()
}
# ------------------------------------------------------------------------------
# 16. Console summary
# ------------------------------------------------------------------------------
cat("\n============================================================\n")
cat("Pure-R exact workflow completed. / R.\n")
cat("Output directory:", normalizePath(output_dir), "\n\n")
print(verification)
cat("\nHigh-yield candidates:\n")
print(actual_high)
cat("\nLow-yield candidates:\n")
print(actual_low)
cat("\nGenerated figures:\n")
cat("- FigF_candidate_metabolite_matrix.svg/png/tiff\n")
cat("- FigG_pathway_module_temporal_response.svg/png/tiff\n")
cat("- FigFG_combined.svg/png/tiff\n")
cat("============================================================\n")
