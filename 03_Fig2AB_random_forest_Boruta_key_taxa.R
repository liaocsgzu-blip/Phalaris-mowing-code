# ============================================================
# Random Forest + Boruta: H vs G key Species screening
# + Boruta:H G
# Input / :species_abundance.csv
# Data format
# Row = Species
# Column = sample relative abundance percentage or proportion
# Example sample names
# RS-H6-1, RT-GCK-1, BS-HCK-5
# Grouping rule
# H group = middle segment starts with H, including H1, H2, ..., HCK
# G group = middle segment starts with G, including G1, G2, ..., GCK
# :RS-H6-1 -> H;RT-GCK-1 -> G
# Main outputs
# 1. All samples H vs G model / H vs G
# 2. Optional separate BS / RS / RT H vs G models
# 3. Boruta importance, MDA ranking, RFCV, model performance
# ============================================================
# ------------------ 0. Clear environment / --------------------
rm(list = ls())
gc()
# ------------------ 1. Load packages / --------------------
packages <- c(
  "randomForest",
  "ggplot2",
  "Boruta",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "data.table"
)
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# ------------------ 2. Parameters / --------------------
set.seed(123)
# Input file
# Keep this R script and species_abundance.csv in the same working directory.
# R species_abundance.csv .
input_file <- "species_abundance.csv"
if (!file.exists(input_file)) {
  candidate_files <- list.files(pattern = "^species_abundance.*\\.csv$", ignore.case = TRUE)
  if (length(candidate_files) > 0) {
    input_file <- candidate_files[1]
    message(" species_abundance.csv,:", input_file)
  } else {
    stop(". species_abundance.csv .")
  }
}
# Output root directory
out_root <- "RF_Boruta_H_vs_G_multi_sample"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
# Feature filtering thresholds
# 0.001 = 0.1% relative abundance after conversion to proportion.
# ,0.001 0.1%.
ra_threshold <- 0.001
detect_threshold <- 0.60
# Random forest and Boruta parameters / Boruta
ntree_num <- 1000
boruta_max_runs <- 300
# Cross-validation parameters
n_repeat <- 10
cv.fold <- 10
step <- 0.9
# Final train/test split
train_fraction <- 0.75
# Whether to additionally run separate models for BS, RS and RT.
# BS,RS,RT .
run_by_sample_type <- TRUE
# ------------------ 3. Helper functions / --------------------
write_csv_utf8 <- function(x, file) {
  utils::write.csv(x, file = file, row.names = FALSE, fileEncoding = "UTF-8")
}
parse_sample_info <- function(sample_names) {
 # Expected format / :SampleType-Treatment-Replicate
 # Examples / :RS-H6-1, RT-GCK-1, BS-HCK-5
  parsed <- stringr::str_match(
    sample_names,
    "^([^-]+)-([HG](?:CK|[0-9]+))-(\\d+)$"
  )
  sample_info <- tibble::tibble(
    Sample = sample_names,
    SampleType = parsed[, 2],
    Treatment = parsed[, 3],
    Replicate = parsed[, 4],
    MainGroup = ifelse(!is.na(parsed[, 3]), substr(parsed[, 3], 1, 1), NA_character_)
  )
  bad_samples <- sample_info$Sample[is.na(sample_info$MainGroup)]
  if (length(bad_samples) > 0) {
    stop(
      paste0(
        " H/G , RS-H6-1,RT-GCK-1,BS-HCK-5 :\n",
        paste(bad_samples, collapse = ", ")
      )
    )
  }
  sample_info <- sample_info %>%
    dplyr::mutate(
      MainGroup = factor(MainGroup, levels = c("G", "H")),
      SampleType = factor(SampleType),
      Treatment = factor(Treatment),
      Replicate = as.integer(Replicate)
    )
  sample_info
}
convert_to_relative_abundance <- function(abund_mat) {
  abund_mat <- as.matrix(abund_mat)
  storage.mode(abund_mat) <- "numeric"
  abund_mat[is.na(abund_mat)] <- 0
  col_sum <- colSums(abund_mat, na.rm = TRUE)
  if (any(col_sum == 0)) {
    stop(
      paste0(
        " 0,:",
        paste(names(col_sum)[col_sum == 0], collapse = ", ")
      )
    )
  }
 # The uploaded file sums to approximately 100 per sample.
 # 100, 100.
  if (all(col_sum > 95 & col_sum < 105)) {
    message(" 100:, 100.")
    rel_mat <- abund_mat / 100
  } else if (all(col_sum > 0.95 & col_sum < 1.05)) {
    message(" 1:,.")
    rel_mat <- abund_mat
  } else {
    message(" 1  100: counts ,.")
    rel_mat <- sweep(abund_mat, 2, col_sum, "/")
  }
  rel_mat[is.na(rel_mat)] <- 0
  as.data.frame(rel_mat, check.names = FALSE)
}
calc_filter_stats <- function(rel_mat, sample_info_use) {
  dplyr::bind_rows(
    lapply(c("G", "H"), function(g) {
      samples_use <- sample_info_use %>%
        dplyr::filter(MainGroup == g) %>%
        dplyr::pull(Sample)
      tibble::tibble(
        FeatureID = rownames(rel_mat),
        MainGroup = g,
        Mean_RA = rowMeans(rel_mat[, samples_use, drop = FALSE], na.rm = TRUE),
        Detection_rate = rowMeans(rel_mat[, samples_use, drop = FALSE] > 0, na.rm = TRUE),
        Pass = Mean_RA >= ra_threshold & Detection_rate >= detect_threshold
      )
    })
  )
}
make_filter_summary <- function(rel_mat, sample_info_use, tax_map) {
  filter_stats_long <- calc_filter_stats(rel_mat, sample_info_use)
  mean_ra_wide <- filter_stats_long %>%
    dplyr::select(FeatureID, MainGroup, Mean_RA) %>%
    tidyr::pivot_wider(
      names_from = MainGroup,
      values_from = Mean_RA,
      names_prefix = "Mean_RA_"
    )
  detect_wide <- filter_stats_long %>%
    dplyr::select(FeatureID, MainGroup, Detection_rate) %>%
    tidyr::pivot_wider(
      names_from = MainGroup,
      values_from = Detection_rate,
      names_prefix = "Detection_rate_"
    )
  pass_wide <- filter_stats_long %>%
    dplyr::select(FeatureID, MainGroup, Pass) %>%
    tidyr::pivot_wider(
      names_from = MainGroup,
      values_from = Pass,
      names_prefix = "Pass_"
    )
  filter_summary <- tax_map %>%
    dplyr::left_join(mean_ra_wide, by = "FeatureID") %>%
    dplyr::left_join(detect_wide, by = "FeatureID") %>%
    dplyr::left_join(pass_wide, by = "FeatureID")
  for (pc in c("Pass_G", "Pass_H")) {
    if (!pc %in% colnames(filter_summary)) {
      filter_summary[[pc]] <- FALSE
    }
    filter_summary[[pc]][is.na(filter_summary[[pc]])] <- FALSE
    filter_summary[[pc]] <- as.logical(filter_summary[[pc]])
  }
  filter_summary <- filter_summary %>%
    dplyr::mutate(
      Pass_any_group = Pass_G | Pass_H,
      Passed_Groups = dplyr::case_when(
        Pass_G & Pass_H ~ "G;H",
        Pass_G & !Pass_H ~ "G",
        !Pass_G & Pass_H ~ "H",
        TRUE ~ ""
      )
    )
  filter_summary
}
safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}
# ------------------ 4. Read data / --------------------
cat(":", input_file, "\n")
raw_data <- data.table::fread(
  input_file,
  data.table = FALSE,
  check.names = FALSE,
  encoding = "UTF-8"
)
colnames(raw_data) <- gsub("^\\ufeff", "", colnames(raw_data))
colnames(raw_data) <- trimws(colnames(raw_data))
cat(":", nrow(raw_data), " ×", ncol(raw_data), "\n")
# ------------------ 5. Identify Species and sample columns / Species --------------------
species_col_candidates <- c("Species", "species", "Taxon", "TaxonLabel")
species_col <- species_col_candidates[species_col_candidates %in% colnames(raw_data)]
if (length(species_col) == 0) {
  species_col <- colnames(raw_data)[1]
  message(" Species , Species:", species_col)
} else {
  species_col <- species_col[1]
  message(" Species :", species_col)
}
raw_data[[species_col]] <- as.character(raw_data[[species_col]])
raw_data[[species_col]][is.na(raw_data[[species_col]]) | raw_data[[species_col]] == ""] <-
  paste0("Unclassified_species_", which(is.na(raw_data[[species_col]]) | raw_data[[species_col]] == ""))
sample_cols <- setdiff(colnames(raw_data), species_col)
sample_cols <- sample_cols[stringr::str_detect(sample_cols, "^[^-]+-[HG](?:CK|[0-9]+)-\\d+$")]
if (length(sample_cols) == 0) {
  stop(". RS-H6-1,RT-GCK-1,BS-HCK-5.")
}
sample_info <- parse_sample_info(sample_cols)
cat("\n:", nrow(sample_info), "\n")
cat("\nH/G :\n")
print(table(sample_info$MainGroup))
cat("\n:\n")
print(table(sample_info$SampleType))
cat("\n × H/G :\n")
print(table(sample_info$SampleType, sample_info$MainGroup))
write_csv_utf8(
  sample_info,
  file.path(out_root, "Sample_group_table_H_vs_G.csv")
)
# ------------------ 6. Prepare abundance matrix / --------------------
abund_raw <- raw_data[, c(species_col, sample_cols), drop = FALSE]
colnames(abund_raw)[colnames(abund_raw) == species_col] <- "Species"
# Convert sample columns to numeric
abund_raw[, sample_cols] <- lapply(abund_raw[, sample_cols, drop = FALSE], function(x) {
  x <- as.numeric(as.character(x))
  x[is.na(x)] <- 0
  x
})
# Merge duplicated Species if present / Species ,
abund_species <- abund_raw %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )
cat("\nSpecies :", nrow(abund_species), "\n")
feature_ids <- paste0("Species_", stringr::str_pad(seq_len(nrow(abund_species)), width = 5, pad = "0"))
tax_map <- tibble::tibble(
  FeatureID = feature_ids,
  Species = abund_species$Species
)
abund_mat <- abund_species %>%
  dplyr::select(dplyr::all_of(sample_cols)) %>%
  as.data.frame(check.names = FALSE)
rownames(abund_mat) <- feature_ids
rel_abund <- convert_to_relative_abundance(abund_mat)
rownames(rel_abund) <- feature_ids
# Output converted relative abundance
rel_out <- tax_map %>%
  dplyr::left_join(
    rel_abund %>%
      tibble::rownames_to_column("FeatureID"),
    by = "FeatureID"
  )
write_csv_utf8(
  rel_out,
  file.path(out_root, "Species_relative_abundance_proportion_all_samples.csv")
)
# ------------------ 7. Core RF + Boruta function / --------------------
run_h_vs_g_rf <- function(samples_use, analysis_name, title_prefix) {
  analysis_name <- safe_filename(analysis_name)
  out_dir <- file.path(out_root, analysis_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  sample_info_use <- sample_info %>%
    dplyr::filter(Sample %in% samples_use) %>%
    dplyr::arrange(match(Sample, samples_use))
  samples_use <- sample_info_use$Sample
  group <- droplevels(sample_info_use$MainGroup)
  cat("\n====================================================\n")
  cat(":", analysis_name, "\n")
  cat(":H vs G\n")
  cat(":", length(samples_use), "\n")
  print(table(group))
  cat("====================================================\n")
  if (length(unique(group)) != 2) {
    stop(" H/G ,:", analysis_name)
  }
  if (min(table(group)) < 3) {
    stop(",:", analysis_name)
  }
 # ---------- 7.1 Feature filtering / ----------
  rel_model <- rel_abund[, samples_use, drop = FALSE]
  filter_summary <- make_filter_summary(rel_model, sample_info_use, tax_map)
  keep_ids <- filter_summary$FeatureID[filter_summary$Pass_any_group]
  cat(" Species :", nrow(rel_model), "\n")
  cat(" Species :", length(keep_ids), "\n")
  cat(
    ": H/G  mean_RA >= ",
    ra_threshold,
    "  detection_rate >= ",
    detect_threshold,
    "\n",
    sep = ""
  )
  write_csv_utf8(
    filter_summary,
    file.path(out_dir, paste0("Species_filter_summary_", analysis_name, ".csv"))
  )
  if (length(keep_ids) < 2) {
    stop(" Species  2, ra_threshold  detect_threshold:", analysis_name)
  }
  rel_model_filtered <- rel_model[keep_ids, , drop = FALSE]
 # ---------- 7.2 Construct RF matrix / ----------
  features <- as.data.frame(
    t(as.matrix(rel_model_filtered)),
    check.names = FALSE
  )
  features <- features[samples_use, , drop = FALSE]
 # Remove zero-variance features / 0
  var_value <- apply(features, 2, var, na.rm = TRUE)
  features <- features[, var_value > 0, drop = FALSE]
  cat(" Species :", ncol(features), "\n")
  if (ncol(features) < 2) {
    stop(" 2,:", analysis_name)
  }
  model_sample_info <- sample_info_use %>%
    dplyr::mutate(ModelGroup = group)
  write_csv_utf8(
    model_sample_info,
    file.path(out_dir, paste0("model_sample_info_", analysis_name, ".csv"))
  )
 # ---------- 7.3 Random forest model / ----------
  set.seed(123)
  rf_model <- randomForest::randomForest(
    x = features,
    y = group,
    importance = TRUE,
    proximity = TRUE,
    ntree = ntree_num
  )
  print(rf_model)
  train_accuracy_all <- mean(predict(rf_model) == group)
 # OOB error
 # English: Out-of-bag error is the internal validation error of Random Forest.
 # :.
  oob_error_all <- tail(rf_model$err.rate[, "OOB"], 1)
  oob_accuracy_all <- 1 - oob_error_all
  class_error_all <- rf_model$confusion[, "class.error"]
  class_error_G_all <- if ("G" %in% names(class_error_all)) unname(class_error_all["G"]) else NA_real_
  class_error_H_all <- if ("H" %in% names(class_error_all)) unname(class_error_all["H"]) else NA_real_
  cat("\nTraining Accuracy using all retained features:", round(train_accuracy_all * 100, 2), "%\n")
  cat("OOB Error using all retained features:", round(oob_error_all, 4), "\n")
  cat("OOB Accuracy using all retained features:", round(oob_accuracy_all, 4), "\n")
  confusion_all <- as.data.frame.matrix(rf_model$confusion)
  confusion_all$TrueGroup <- rownames(confusion_all)
  confusion_all <- confusion_all %>%
    dplyr::select(TrueGroup, dplyr::everything())
  write_csv_utf8(
    confusion_all,
    file.path(out_dir, paste0("RF_OOB_confusion_matrix_all_retained_features_", analysis_name, ".csv"))
  )
  capture.output(
    rf_model,
    file = file.path(out_dir, paste0("random_forest_model_", analysis_name, ".txt"))
  )
 # ---------- 7.4 Boruta / Boruta ----------
  set.seed(123)
  boruta_result <- Boruta::Boruta(
    x = features,
    y = group,
    doTrace = 1,
    maxRuns = boruta_max_runs
  )
  pdf(
    file = file.path(out_dir, paste0("01_Boruta_Variable_Importance_", analysis_name, ".pdf")),
    width = 12,
    height = 7
  )
  plot(
    boruta_result,
    las = 2,
    cex.axis = 0.6,
    main = paste0(title_prefix, ": Boruta Variable Importance")
  )
  dev.off()
  png(
    file = file.path(out_dir, paste0("01_Boruta_Variable_Importance_", analysis_name, ".png")),
    width = 3600,
    height = 2100,
    res = 300
  )
  plot(
    boruta_result,
    las = 2,
    cex.axis = 0.6,
    main = paste0(title_prefix, ": Boruta Variable Importance")
  )
  dev.off()
  boruta_stats <- Boruta::attStats(boruta_result)
  boruta_stats$FeatureID <- rownames(boruta_stats)
  boruta_stats$decision <- as.character(boruta_stats$decision)
  write_csv_utf8(
    boruta_stats,
    file.path(out_dir, paste0("boruta_statistics_", analysis_name, ".csv"))
  )
 # ---------- 7.5 MDA + mean abundance + fold change / MDA, ----------
  importance_data <- randomForest::importance(rf_model, type = 1)
  mda_df <- data.frame(
    FeatureID = rownames(importance_data),
    MeanDecreaseAccuracy = importance_data[, "MeanDecreaseAccuracy"],
    stringsAsFactors = FALSE
  )
  mean_g <- colMeans(features[group == "G", , drop = FALSE], na.rm = TRUE)
  mean_h <- colMeans(features[group == "H", , drop = FALSE], na.rm = TRUE)
  pseudo <- 1e-9
  mean_df <- data.frame(
    FeatureID = names(mean_g),
    Mean_G = as.numeric(mean_g),
    Mean_H = as.numeric(mean_h),
    stringsAsFactors = FALSE
  )
  mean_df$FoldChange_H_vs_G <- (mean_df$Mean_H + pseudo) / (mean_df$Mean_G + pseudo)
  mean_df$log2FC_H_vs_G <- log2(mean_df$FoldChange_H_vs_G)
  mean_df$ChangeDirection <- ifelse(
    mean_df$Mean_H > mean_df$Mean_G,
    "Higher in H",
    ifelse(mean_df$Mean_H < mean_df$Mean_G, "Lower in H", "No change")
  )
  mda_df <- mda_df %>%
    dplyr::left_join(
      boruta_stats[, c("FeatureID", "decision")],
      by = "FeatureID"
    ) %>%
    dplyr::left_join(mean_df, by = "FeatureID") %>%
    dplyr::left_join(tax_map, by = "FeatureID")
  mda_df$decision <- as.character(mda_df$decision)
  mda_df$decision[is.na(mda_df$decision)] <- "Not_tested"
  mda_df$decision <- factor(
    mda_df$decision,
    levels = c("Confirmed", "Tentative", "Rejected", "Not_tested")
  )
  mda_df$significance <- ifelse(mda_df$decision == "Confirmed", "*", "")
  mda_df <- mda_df %>%
    dplyr::arrange(dplyr::desc(MeanDecreaseAccuracy))
  write_csv_utf8(
    mda_df,
    file.path(out_dir, paste0("all_mda_species_with_boruta_", analysis_name, ".csv"))
  )
  top_mda <- head(mda_df, 20)
  top_mda$PlotLabel <- stringr::str_wrap(top_mda$Species, width = 45)
  top_mda$PlotLabel[is.na(top_mda$PlotLabel) | top_mda$PlotLabel == ""] <-
    top_mda$FeatureID[is.na(top_mda$PlotLabel) | top_mda$PlotLabel == ""]
  top_mda$PlotLabel <- make.unique(top_mda$PlotLabel)
  write_csv_utf8(
    top_mda,
    file.path(out_dir, paste0("top20_mda_species_with_boruta_", analysis_name, ".csv"))
  )
 # ---------- 7.6 Top20 MDA plot / Top20 MDA ----------
  p_mda <- ggplot2::ggplot(
    top_mda,
    ggplot2::aes(
      x = reorder(PlotLabel, MeanDecreaseAccuracy),
      y = MeanDecreaseAccuracy,
      fill = decision
    )
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(
      ggplot2::aes(label = significance),
      hjust = -0.3,
      size = 5
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Species",
      y = "Mean Decrease Accuracy",
      title = paste0(title_prefix, ": Top 20 species by MDA"),
      fill = "Boruta decision"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(
      values = c(
        "Confirmed" = "seagreen",
        "Tentative" = "gold",
        "Rejected" = "gray",
        "Not_tested" = "gray80"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.text.x = ggplot2::element_text(size = 10),
      legend.position = "right"
    )
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("02_Top20_MDA_", analysis_name, ".pdf")),
    plot = p_mda,
    width = 10,
    height = 7
  )
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("02_Top20_MDA_", analysis_name, ".png")),
    plot = p_mda,
    width = 10,
    height = 7,
    dpi = 300
  )
 # ---------- 7.7 Repeated RFCV / ----------
  actual_cv_fold <- min(cv.fold, as.integer(min(table(group))))
  cat("\n:", actual_cv_fold, "\n")
  cv_all <- list()
  set.seed(123)
  result <- randomForest::rfcv(
    trainx = features,
    trainy = group,
    cv.fold = actual_cv_fold,
    step = step
  )
  cv_all[[1]] <- result$error.cv
  nvar_seq <- result$n.var
  for (i in 2:n_repeat) {
    set.seed(123 + i)
    result_i <- randomForest::rfcv(
      trainx = features,
      trainy = group,
      cv.fold = actual_cv_fold,
      step = step
    )
    cv_all[[i]] <- result_i$error.cv
  }
  cv_df <- data.frame(NumVariables = nvar_seq)
  for (i in 1:n_repeat) {
    cv_df[[paste0("error.", i)]] <- cv_all[[i]]
  }
  error_cols <- paste0("error.", 1:n_repeat)
  cv_df$ErrorMean <- rowMeans(cv_df[, error_cols, drop = FALSE], na.rm = TRUE)
  optimal_n <- cv_df$NumVariables[which.min(cv_df$ErrorMean)]
  write_csv_utf8(
    cv_df,
    file.path(out_dir, paste0("rfcv_error_summary_", analysis_name, ".csv"))
  )
  cv_long <- cv_df %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(error_cols),
      names_to = "Replicate",
      values_to = "Error"
    )
  p_cv <- ggplot2::ggplot(
    cv_long,
    ggplot2::aes(
      x = NumVariables,
      y = Error,
      group = Replicate
    )
  ) +
    ggplot2::geom_line(color = "gray", linewidth = 0.5) +
    ggplot2::geom_line(
      data = cv_df,
      ggplot2::aes(x = NumVariables, y = ErrorMean),
      inherit.aes = FALSE,
      linewidth = 1.2,
      color = "black"
    ) +
    ggplot2::geom_vline(
      xintercept = optimal_n,
      linetype = "dashed",
      color = "red"
    ) +
    ggplot2::annotate(
      "text",
      x = optimal_n,
      y = max(cv_df$ErrorMean, na.rm = TRUE),
      label = paste("Optimal =", optimal_n),
      vjust = -1
    ) +
    ggplot2::labs(
      title = paste0(title_prefix, ": ", actual_cv_fold, "-fold cross-validation repeated ", n_repeat, " times"),
      x = "Number of Variables",
      y = "Cross-validation Error"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      axis.text = ggplot2::element_text(size = 10)
    )
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("03_RFCV_", analysis_name, ".pdf")),
    plot = p_cv,
    width = 8,
    height = 6
  )
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("03_RFCV_", analysis_name, ".png")),
    plot = p_cv,
    width = 8,
    height = 6,
    dpi = 300
  )
 # ---------- 7.8 Final train/test model / - ----------
  important_vars <- mda_df$FeatureID[order(mda_df$MeanDecreaseAccuracy, decreasing = TRUE)]
  important_vars <- important_vars[1:min(optimal_n, length(important_vars))]
  reduced_features <- features[, important_vars, drop = FALSE]
 # OOB performance of the reduced feature model fitted on all samples.
 # English: This OOB metric is a stability reference after RFCV feature selection.
 # : RFCV ,,.
 # :,.
  set.seed(123)
  rf_reduced_all_model <- randomForest::randomForest(
    x = reduced_features,
    y = group,
    importance = TRUE,
    proximity = TRUE,
    ntree = ntree_num
  )
  oob_error_reduced_all <- tail(rf_reduced_all_model$err.rate[, "OOB"], 1)
  oob_accuracy_reduced_all <- 1 - oob_error_reduced_all
  class_error_reduced_all <- rf_reduced_all_model$confusion[, "class.error"]
  class_error_G_reduced_all <- if ("G" %in% names(class_error_reduced_all)) unname(class_error_reduced_all["G"]) else NA_real_
  class_error_H_reduced_all <- if ("H" %in% names(class_error_reduced_all)) unname(class_error_reduced_all["H"]) else NA_real_
  confusion_reduced_all <- as.data.frame.matrix(rf_reduced_all_model$confusion)
  confusion_reduced_all$TrueGroup <- rownames(confusion_reduced_all)
  confusion_reduced_all <- confusion_reduced_all %>%
    dplyr::select(TrueGroup, dplyr::everything())
  write_csv_utf8(
    confusion_reduced_all,
    file.path(out_dir, paste0("RF_OOB_confusion_matrix_reduced_features_", analysis_name, ".csv"))
  )
  capture.output(
    rf_reduced_all_model,
    file = file.path(out_dir, paste0("random_forest_reduced_features_model_", analysis_name, ".txt"))
  )
  set.seed(123)
  train_index <- unlist(
    tapply(seq_along(group), group, function(idx) {
      sample(idx, size = max(1, floor(length(idx) * train_fraction)))
    })
  )
  train_data <- reduced_features[train_index, , drop = FALSE]
  train_label <- group[train_index]
  test_data <- reduced_features[-train_index, , drop = FALSE]
  test_label <- group[-train_index]
  set.seed(123)
  rf_test_model <- randomForest::randomForest(
    x = train_data,
    y = train_label,
    ntree = ntree_num
  )
  train_pred <- predict(rf_test_model, train_data)
  test_pred <- predict(rf_test_model, test_data)
  acc_train <- mean(train_pred == train_label)
  acc_test <- mean(test_pred == test_label)
  err_train <- 1 - acc_train
  err_test <- 1 - acc_test
  feat_num <- ncol(reduced_features)
  cat("\n==== :", analysis_name, " ====\n")
  cat("Accuracy Train Set:", round(acc_train, 4), "\n")
  cat("Accuracy Test Set:", round(acc_test, 4), "\n")
  cat("Feature Number:", feat_num, "\n")
  cat("Error Rate Train Set:", round(err_train, 4), "\n")
  cat("Error Rate Test Set:", round(err_test, 4), "\n")
  cat("=========================================\n")
  metrics <- data.frame(
    Analysis = analysis_name,
    Comparison = "H_vs_G",
    TotalSampleNumber = nrow(features),
    G_SampleNumber = sum(group == "G"),
    H_SampleNumber = sum(group == "H"),
    SpeciesNumber_AfterFilter = ncol(features),
    TrainingAccuracy_AllRetainedFeatures = round(train_accuracy_all, 4),
    OOB_Error_AllRetainedFeatures = round(oob_error_all, 4),
    OOB_Accuracy_AllRetainedFeatures = round(oob_accuracy_all, 4),
    OOB_ClassError_G_AllRetainedFeatures = round(class_error_G_all, 4),
    OOB_ClassError_H_AllRetainedFeatures = round(class_error_H_all, 4),
    OOB_Error_ReducedFeatures_AllSamples = round(oob_error_reduced_all, 4),
    OOB_Accuracy_ReducedFeatures_AllSamples = round(oob_accuracy_reduced_all, 4),
    OOB_ClassError_G_ReducedFeatures = round(class_error_G_reduced_all, 4),
    OOB_ClassError_H_ReducedFeatures = round(class_error_H_reduced_all, 4),
    Accuracy_TrainSet_FinalModel = round(acc_train, 4),
    Accuracy_TestSet_FinalModel = round(acc_test, 4),
    Feature_Number_FinalModel = feat_num,
    Error_Rate_TrainSet_FinalModel = round(err_train, 4),
    Error_Rate_TestSet_FinalModel = round(err_test, 4),
    Actual_CV_Fold = actual_cv_fold,
    Optimal_Feature_Number_RFCV = optimal_n,
    stringsAsFactors = FALSE
  )
  write_csv_utf8(
    metrics,
    file.path(out_dir, paste0("RF_model_performance_summary_", analysis_name, ".csv"))
  )
  final_selected_species <- mda_df[mda_df$FeatureID %in% important_vars, ]
  final_selected_species <- final_selected_species[
    match(important_vars, final_selected_species$FeatureID),
  ]
  write_csv_utf8(
    final_selected_species,
    file.path(out_dir, paste0("final_selected_species_", analysis_name, ".csv"))
  )
  cat("\n:", analysis_name, "\n")
  cat(":", out_dir, "\n")
  list(
    analysis_name = analysis_name,
    top_mda = top_mda,
    all_mda = mda_df,
    metrics = metrics,
    final_selected_species = final_selected_species,
    cv_df = cv_df
  )
}
# ------------------ 8. Run models / --------------------
results <- list()
# 8.1 All samples / H vs G
results[["All_samples_H_vs_G"]] <- run_h_vs_g_rf(
  samples_use = sample_info$Sample,
  analysis_name = "All_samples_H_vs_G",
  title_prefix = "All samples H-responsive vs G taxa"
)
# 8.2 Optional separate models by sample type
if (isTRUE(run_by_sample_type)) {
  for (st in levels(sample_info$SampleType)) {
    st_samples <- sample_info %>%
      dplyr::filter(SampleType == st) %>%
      dplyr::pull(Sample)
    results[[paste0(st, "_H_vs_G")]] <- run_h_vs_g_rf(
      samples_use = st_samples,
      analysis_name = paste0(st, "_H_vs_G"),
      title_prefix = paste0(st, " samples H-responsive vs G taxa")
    )
  }
}
# ------------------ 9. Combined outputs / --------------------
metrics_combined <- dplyr::bind_rows(lapply(results, function(x) x$metrics))
write_csv_utf8(
  metrics_combined,
  file.path(out_root, "RF_model_performance_summary_combined_H_vs_G.csv")
)
top20_combined <- dplyr::bind_rows(lapply(results, function(x) {
  x$top_mda %>%
    dplyr::mutate(Analysis = x$analysis_name)
}))
write_csv_utf8(
  top20_combined,
  file.path(out_root, "Top20_MDA_combined_H_vs_G.csv")
)
final_selected_combined <- dplyr::bind_rows(lapply(results, function(x) {
  x$final_selected_species %>%
    dplyr::mutate(Analysis = x$analysis_name)
}))
write_csv_utf8(
  final_selected_combined,
  file.path(out_root, "Final_selected_species_combined_H_vs_G.csv")
)
# ------------------ 10. Finish / --------------------
cat("\n====================================================\n")
cat(".\n\n")
cat(":", input_file, "\n")
cat(":H vs G\n")
cat("H/G :, RS-H6-1  H,RT-GCK-1  G.\n")
cat(":, counts,.\n")
cat(
  ": H/G  mean_RA >= ",
  ra_threshold,
  "  detection_rate >= ",
  detect_threshold,
  "\n",
  sep = ""
)
cat(":", out_root, "\n")
cat("====================================================\n")
