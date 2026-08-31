# ============================================================
# Key Species Profile Contribution Rate
# BS / RS / RT
# 1) species_abundance.csv species_abundance.csv
# 2) final_selected_species_All_samples_H_vs_G.csv
# All samples 43 Species ,
# BS,RS,RT H vs G .
# = sum(|Mean_H - Mean_G|)
# MDA = sum(|Mean_H - Mean_G| * MeanDecreaseAccuracy)
# .
# .
# "".
# ============================================================
rm(list = ls())
gc()
# ------------------ 1. --------------------
packages <- c("data.table", "dplyr", "tidyr", "stringr", "ggplot2", "tibble")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
set.seed(123)
# ------------------ 2. --------------------
# ,.
species_file <- "species_abundance.csv"
key_file <- file.path(
  "RF_Boruta_H_vs_G_multi_sample",
  "All_samples_H_vs_G",
  "final_selected_species_All_samples_H_vs_G.csv"
)
# (2) (1),
if (!file.exists(species_file)) {
  species_candidates <- list.files(pattern = "^Species.*\\.csv$", ignore.case = FALSE)
  if (length(species_candidates) > 0) {
    species_file <- species_candidates[1]
    message(" Species :", species_file)
  } else {
    stop(" species_abundance.csv / species_abundance.csv.")
  }
}
if (!file.exists(key_file)) {
  key_candidates <- list.files(
    pattern = "^final_selected_species_All_samples_H_vs_G.*\\.csv$",
    ignore.case = FALSE,
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(key_candidates) > 0) {
    key_file <- key_candidates[1]
    message(":", key_file)
  } else {
    stop(" final_selected_species_All_samples_H_vs_G.csv / final_selected_species_All_samples_H_vs_G(1).csv.")
  }
}
out_dir <- "Key_species_profile_contribution_rate_ONLY_two_files"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv_utf8 <- function(x, file) {
  utils::write.csv(x, file = file, row.names = FALSE, fileEncoding = "UTF-8")
}
# ------------------ 3. --------------------
cat(":", species_file, "\n")
species_raw <- data.table::fread(
  species_file,
  data.table = FALSE,
  check.names = FALSE,
  encoding = "UTF-8"
)
cat(":", key_file, "\n")
key_species <- data.table::fread(
  key_file,
  data.table = FALSE,
  check.names = FALSE,
  encoding = "UTF-8"
)
colnames(species_raw) <- trimws(gsub("^\\ufeff", "", colnames(species_raw)))
colnames(key_species) <- trimws(gsub("^\\ufeff", "", colnames(key_species)))
if (!"Species" %in% colnames(species_raw)) {
  stop("Species  Species .")
}
if (!"Species" %in% colnames(key_species)) {
  stop(" Species .")
}
if (!"MeanDecreaseAccuracy" %in% colnames(key_species)) {
  stop(" MeanDecreaseAccuracy , MDA .")
}
if (!"ChangeDirection" %in% colnames(key_species)) {
  stop(" ChangeDirection , H  G .")
}
# Species,
key_species <- key_species %>%
  dplyr::mutate(
    Species = as.character(Species),
    MeanDecreaseAccuracy = as.numeric(MeanDecreaseAccuracy),
    ChangeDirection = as.character(ChangeDirection)
  ) %>%
  dplyr::distinct(Species, .keep_all = TRUE)
cat(":", nrow(key_species), "\n")
# ------------------ 4. --------------------
sample_cols <- setdiff(colnames(species_raw), "Species")
sample_cols <- sample_cols[
  stringr::str_detect(sample_cols, "^[^-]+-[HG](?:CK|[0-9]+)-[0-9]+$")
]
if (length(sample_cols) == 0) {
  stop(". BS-H6-1,RS-GCK-1,RT-H16-5.")
}
sample_info <- tibble::tibble(Sample = sample_cols) %>%
  dplyr::mutate(
    Profile = stringr::str_match(Sample, "^([^-]+)-([HG](?:CK|[0-9]+))-([0-9]+)$")[, 2],
    Treatment = stringr::str_match(Sample, "^([^-]+)-([HG](?:CK|[0-9]+))-([0-9]+)$")[, 3],
    Replicate = as.integer(stringr::str_match(Sample, "^([^-]+)-([HG](?:CK|[0-9]+))-([0-9]+)$")[, 4]),
    Group = substr(Treatment, 1, 1)
  ) %>%
  dplyr::filter(!is.na(Profile), Group %in% c("G", "H"))
if (nrow(sample_info) == 0) {
  stop("..")
}
# ;,
profile_order <- c("BS", "RS", "RT")
profile_levels <- unique(c(profile_order[profile_order %in% unique(sample_info$Profile)],
                           sort(setdiff(unique(sample_info$Profile), profile_order))))
sample_info$Profile <- factor(sample_info$Profile, levels = profile_levels)
sample_info$Group <- factor(sample_info$Group, levels = c("G", "H"))
cat("\n:", nrow(sample_info), "\n")
cat("\n × H/G :\n")
print(table(sample_info$Profile, sample_info$Group))
write_csv_utf8(sample_info, file.path(out_dir, "Sample_info_parsed.csv"))
# ------------------ 5. --------------------
species_raw[, sample_info$Sample] <- lapply(species_raw[, sample_info$Sample, drop = FALSE], function(x) {
  x <- as.numeric(as.character(x))
  x[is.na(x)] <- 0
  x
})
# Species ,
species_abund <- species_raw %>%
  dplyr::select(Species, dplyr::all_of(sample_info$Sample)) %>%
  dplyr::mutate(Species = as.character(Species)) %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(sample_info$Sample), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )
abund_mat <- species_abund %>%
  tibble::column_to_rownames("Species") %>%
  as.data.frame(check.names = FALSE)
sample_sum <- colSums(abund_mat[, sample_info$Sample, drop = FALSE], na.rm = TRUE)
if (all(sample_sum > 95 & sample_sum < 105)) {
  cat("\n 100,, 100.\n")
  rel_mat <- abund_mat / 100
} else if (all(sample_sum > 0.95 & sample_sum < 1.05)) {
  cat("\n 1,.\n")
  rel_mat <- abund_mat
} else {
  cat("\n 100  1,.\n")
  rel_mat <- sweep(abund_mat, 2, sample_sum, "/")
}
rel_mat[is.na(rel_mat)] <- 0
missing_species <- setdiff(key_species$Species, rownames(rel_mat))
if (length(missing_species) > 0) {
  write_csv_utf8(
    data.frame(MissingSpecies = missing_species),
    file.path(out_dir, "Missing_key_species_in_species_table.csv")
  )
  stop(
    paste0(
      " ", length(missing_species),
      "  Species . Missing_key_species_in_species_table.csv."
    )
  )
}
rel_key <- rel_mat[key_species$Species, sample_info$Sample, drop = FALSE]
# ------------------ 6. --------------------
calc_profile_contribution <- function(sp, profile_value) {
  h_samples <- sample_info %>%
    dplyr::filter(Profile == profile_value, Group == "H") %>%
    dplyr::pull(Sample)
  g_samples <- sample_info %>%
    dplyr::filter(Profile == profile_value, Group == "G") %>%
    dplyr::pull(Sample)
  if (length(h_samples) == 0 || length(g_samples) == 0) {
    return(NULL)
  }
  mean_h <- mean(as.numeric(rel_key[sp, h_samples]), na.rm = TRUE)
  mean_g <- mean(as.numeric(rel_key[sp, g_samples]), na.rm = TRUE)
  delta <- mean_h - mean_g
  abs_delta <- abs(delta)
  key_row <- key_species[key_species$Species == sp, , drop = FALSE]
  data.frame(
    Species = sp,
    Profile = as.character(profile_value),
    Mean_G = mean_g,
    Mean_H = mean_h,
    Delta_H_minus_G = delta,
    AbsDelta = abs_delta,
    MeanDecreaseAccuracy = key_row$MeanDecreaseAccuracy[1],
    MDA_weighted_AbsDelta = abs_delta * key_row$MeanDecreaseAccuracy[1],
    GlobalChangeDirection = key_row$ChangeDirection[1],
    stringsAsFactors = FALSE
  )
}
contribution_detail <- dplyr::bind_rows(
  lapply(key_species$Species, function(sp) {
    dplyr::bind_rows(
      lapply(levels(sample_info$Profile), function(pf) {
        calc_profile_contribution(sp, pf)
      })
    )
  })
)
contribution_detail <- contribution_detail %>%
  dplyr::mutate(
    DirectionGroup = dplyr::case_when(
      GlobalChangeDirection == "Higher in H" ~ "H_enriched_key_taxa",
      GlobalChangeDirection == "Lower in H" ~ "G_enriched_key_taxa",
      TRUE ~ "Other"
    ),
    Profile = factor(Profile, levels = profile_levels)
  ) %>%
  dplyr::arrange(Species, Profile)
write_csv_utf8(
  contribution_detail,
  file.path(out_dir, "Key_species_contribution_by_taxon_and_profile.csv")
)
# ------------------ 7. --------------------
profile_contribution <- contribution_detail %>%
  dplyr::group_by(Profile) %>%
  dplyr::summarise(
    KeySpeciesNumber = dplyr::n_distinct(Species),
    UnweightedContribution = sum(AbsDelta, na.rm = TRUE),
    MDA_weighted_Contribution = sum(MDA_weighted_AbsDelta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    UnweightedContributionRate_percent =
      UnweightedContribution / sum(UnweightedContribution, na.rm = TRUE) * 100,
    MDA_weighted_ContributionRate_percent =
      MDA_weighted_Contribution / sum(MDA_weighted_Contribution, na.rm = TRUE) * 100
  ) %>%
  dplyr::arrange(dplyr::desc(MDA_weighted_ContributionRate_percent))
write_csv_utf8(
  profile_contribution,
  file.path(out_dir, "Profile_contribution_rate_summary.csv")
)
# MDA ,
write_csv_utf8(
  profile_contribution %>%
    dplyr::select(Profile, KeySpeciesNumber, UnweightedContribution, UnweightedContributionRate_percent),
  file.path(out_dir, "Profile_contribution_rate_unweighted.csv")
)
write_csv_utf8(
  profile_contribution %>%
    dplyr::select(Profile, KeySpeciesNumber, MDA_weighted_Contribution, MDA_weighted_ContributionRate_percent),
  file.path(out_dir, "Profile_contribution_rate_MDA_weighted.csv")
)
# ------------------ 8. H / G --------------------
directional_contribution <- contribution_detail %>%
  dplyr::group_by(DirectionGroup, Profile) %>%
  dplyr::summarise(
    KeySpeciesNumber = dplyr::n_distinct(Species),
    UnweightedContribution = sum(AbsDelta, na.rm = TRUE),
    MDA_weighted_Contribution = sum(MDA_weighted_AbsDelta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(DirectionGroup) %>%
  dplyr::mutate(
    UnweightedContributionRate_within_direction_percent =
      UnweightedContribution / sum(UnweightedContribution, na.rm = TRUE) * 100,
    MDA_weighted_ContributionRate_within_direction_percent =
      MDA_weighted_Contribution / sum(MDA_weighted_Contribution, na.rm = TRUE) * 100
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(DirectionGroup, Profile)
write_csv_utf8(
  directional_contribution,
  file.path(out_dir, "Directional_profile_contribution_rate.csv")
)
# ------------------ 9. --------------------
taxon_profile_share <- contribution_detail %>%
  dplyr::group_by(Species) %>%
  dplyr::mutate(
    TotalAbsDelta_across_profiles = sum(AbsDelta, na.rm = TRUE),
    ProfileShare_percent = ifelse(
      TotalAbsDelta_across_profiles > 0,
      AbsDelta / TotalAbsDelta_across_profiles * 100,
      NA_real_
    )
  ) %>%
  dplyr::ungroup()
dominant_profile_by_taxon <- taxon_profile_share %>%
  dplyr::group_by(Species) %>%
  dplyr::slice_max(order_by = AbsDelta, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    Species,
    DominantProfile = Profile,
    DominantProfileShare_percent = ProfileShare_percent,
    GlobalChangeDirection,
    MeanDecreaseAccuracy,
    Mean_G,
    Mean_H,
    Delta_H_minus_G,
    AbsDelta
  ) %>%
  dplyr::arrange(DominantProfile, dplyr::desc(DominantProfileShare_percent))
write_csv_utf8(
  taxon_profile_share,
  file.path(out_dir, "Key_species_profile_share_by_taxon.csv")
)
write_csv_utf8(
  dominant_profile_by_taxon,
  file.path(out_dir, "Dominant_profile_by_key_species.csv")
)
dominant_profile_count <- dominant_profile_by_taxon %>%
  dplyr::count(DominantProfile, name = "DominantKeySpeciesNumber") %>%
  dplyr::mutate(
    DominantKeySpeciesRate_percent =
      DominantKeySpeciesNumber / sum(DominantKeySpeciesNumber) * 100
  ) %>%
  dplyr::arrange(dplyr::desc(DominantKeySpeciesNumber))
write_csv_utf8(
  dominant_profile_count,
  file.path(out_dir, "Dominant_profile_taxa_count.csv")
)
# ------------------ 10. --------------------
top_weighted_profile <- profile_contribution$Profile[1]
top_weighted_rate <- profile_contribution$MDA_weighted_ContributionRate_percent[1]
top_unweighted_profile <- profile_contribution$Profile[
  which.max(profile_contribution$UnweightedContributionRate_percent)
]
top_unweighted_rate <- max(profile_contribution$UnweightedContributionRate_percent)
interpretation_summary <- data.frame(
  Metric = c(
    "Key_species_number",
    "Top_profile_by_MDA_weighted_contribution",
    "Top_profile_MDA_weighted_contribution_rate_percent",
    "Top_profile_by_unweighted_contribution",
    "Top_profile_unweighted_contribution_rate_percent"
  ),
  Value = c(
    nrow(key_species),
    as.character(top_weighted_profile),
    round(top_weighted_rate, 3),
    as.character(top_unweighted_profile),
    round(top_unweighted_rate, 3)
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(
  interpretation_summary,
  file.path(out_dir, "Interpretation_summary_profile_contribution.csv")
)
# ------------------ 11. --------------------
# 11.1 MDA
p_weighted <- ggplot2::ggplot(
  profile_contribution,
  ggplot2::aes(
    x = reorder(as.character(Profile), -MDA_weighted_ContributionRate_percent),
    y = MDA_weighted_ContributionRate_percent
  )
) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(round(MDA_weighted_ContributionRate_percent, 1), "%")),
    vjust = -0.4,
    size = 4
  ) +
  ggplot2::labs(
    x = "Soil profile",
    y = "MDA-weighted contribution rate (%)",
    title = "Profile contribution of all-sample key species"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5)
  ) +
  ggplot2::ylim(0, max(profile_contribution$MDA_weighted_ContributionRate_percent) * 1.2)
ggplot2::ggsave(
  file.path(out_dir, "01_Profile_contribution_rate_MDA_weighted.pdf"),
  p_weighted,
  width = 5.5,
  height = 4.5
)
ggplot2::ggsave(
  file.path(out_dir, "01_Profile_contribution_rate_MDA_weighted.png"),
  p_weighted,
  width = 5.5,
  height = 4.5,
  dpi = 300
)
# 11.2
p_unweighted <- ggplot2::ggplot(
  profile_contribution,
  ggplot2::aes(
    x = reorder(as.character(Profile), -UnweightedContributionRate_percent),
    y = UnweightedContributionRate_percent
  )
) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(round(UnweightedContributionRate_percent, 1), "%")),
    vjust = -0.4,
    size = 4
  ) +
  ggplot2::labs(
    x = "Soil profile",
    y = "Unweighted contribution rate (%)",
    title = "Profile contribution based on absolute H-G differences"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5)
  ) +
  ggplot2::ylim(0, max(profile_contribution$UnweightedContributionRate_percent) * 1.2)
ggplot2::ggsave(
  file.path(out_dir, "02_Profile_contribution_rate_unweighted.pdf"),
  p_unweighted,
  width = 5.5,
  height = 4.5
)
ggplot2::ggsave(
  file.path(out_dir, "02_Profile_contribution_rate_unweighted.png"),
  p_unweighted,
  width = 5.5,
  height = 4.5,
  dpi = 300
)
# 11.3
p_direction <- ggplot2::ggplot(
  directional_contribution,
  ggplot2::aes(
    x = Profile,
    y = MDA_weighted_ContributionRate_within_direction_percent,
    fill = DirectionGroup
  )
) +
  ggplot2::geom_col(position = "dodge", width = 0.65) +
  ggplot2::labs(
    x = "Soil profile",
    y = "MDA-weighted contribution rate within direction (%)",
    title = "Directional profile contribution of key species",
    fill = "Key-taxa direction"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5)
  )
ggplot2::ggsave(
  file.path(out_dir, "03_Directional_profile_contribution_rate.pdf"),
  p_direction,
  width = 7,
  height = 4.8
)
ggplot2::ggsave(
  file.path(out_dir, "03_Directional_profile_contribution_rate.png"),
  p_direction,
  width = 7,
  height = 4.8,
  dpi = 300
)
# 11.4 H-G
heatmap_df <- contribution_detail %>%
  dplyr::mutate(
    Species = factor(Species, levels = key_species$Species),
    Profile = factor(Profile, levels = profile_levels)
  )
p_heatmap <- ggplot2::ggplot(
  heatmap_df,
  ggplot2::aes(
    x = Profile,
    y = Species,
    fill = Delta_H_minus_G
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  ggplot2::labs(
    x = "Soil profile",
    y = "Key species",
    fill = "H - G",
    title = "Profile-specific H-G differences of key species"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    axis.text.y = ggplot2::element_text(size = 6)
  )
ggplot2::ggsave(
  file.path(out_dir, "04_H_minus_G_delta_heatmap_key_species.pdf"),
  p_heatmap,
  width = 6.5,
  height = 10
)
ggplot2::ggsave(
  file.path(out_dir, "04_H_minus_G_delta_heatmap_key_species.png"),
  p_heatmap,
  width = 6.5,
  height = 10,
  dpi = 300
)
# 11.5
p_count <- ggplot2::ggplot(
  dominant_profile_count,
  ggplot2::aes(
    x = DominantProfile,
    y = DominantKeySpeciesNumber
  )
) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = DominantKeySpeciesNumber),
    vjust = -0.4,
    size = 4
  ) +
  ggplot2::labs(
    x = "Dominant soil profile",
    y = "Number of key species",
    title = "Dominant profile of individual key species"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5)
  ) +
  ggplot2::ylim(0, max(dominant_profile_count$DominantKeySpeciesNumber) * 1.2)
ggplot2::ggsave(
  file.path(out_dir, "05_Dominant_profile_taxa_count.pdf"),
  p_count,
  width = 5.5,
  height = 4.5
)
ggplot2::ggsave(
  file.path(out_dir, "05_Dominant_profile_taxa_count.png"),
  p_count,
  width = 5.5,
  height = 4.5,
  dpi = 300
)
# ------------------ 12. --------------------
cat("\n====================================================\n")
cat(".\n\n")
cat(":\n")
cat("1) ", species_file, "\n", sep = "")
cat("2) ", key_file, "\n\n", sep = "")
cat(":", nrow(key_species), "\n")
cat(":", paste(levels(sample_info$Profile), collapse = ", "), "\n\n")
cat("MDA :\n")
print(
  profile_contribution %>%
    dplyr::select(Profile, MDA_weighted_ContributionRate_percent)
)
cat("\n:\n")
print(
  profile_contribution %>%
    dplyr::select(Profile, UnweightedContributionRate_percent)
)
cat("\nMDA :", as.character(top_weighted_profile),
    ", = ", round(top_weighted_rate, 2), "%\n", sep = "")
cat("\n:", out_dir, "\n")
cat("====================================================\n")
