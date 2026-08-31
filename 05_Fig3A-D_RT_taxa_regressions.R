# -*- coding: UTF-8 -*-
# ============================================================
# 4RTRS
# 1. main_plant_soil_data.xlsx
# 2. species_abundance.csv
# 4 × 4 = 16
# Yield
# RS_TN
# RS_TK
# RS_AP
# Response = β0 + β1 × Relative abundance + ε
# GH
# .
# 1. 16PDF
# 2. 16SVG
# Export publication-ready PDF files
# Export publication-ready SVG files
# 5.
# 6.
# ,PFDR.
# ============================================================
# ============================================================
# 0. R
# ============================================================
required_packages <- c(
  "readxl",
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "tibble",
  "ggplot2",
  "patchwork"
)
missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "R:\n\n",
      "install.packages(c(",
      paste(
        sprintf(
          '"%s"',
          missing_packages
        ),
        collapse = ", "
      ),
      "))"
    )
  )
}
suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
# ============================================================
# 1.
# ============================================================
work_dir <- "."
excel_file <- file.path(
  work_dir,
  "main_plant_soil_data.xlsx"
)
species_file <- file.path(
  work_dir,
  "species_abundance.csv"
)
if (!file.exists(excel_file)) {
  stop(
    paste0(
      "Excel:\n",
      excel_file
    )
  )
}
if (!file.exists(species_file)) {
  stop(
    paste0(
      ":\n",
      species_file
    )
  )
}
# ============================================================
# 2.
# ============================================================
output_dir <- file.path(
  work_dir,
  "Four_RT_taxa_Yield_RS_nutrients_regression"
)
individual_pdf_dir <- file.path(
  output_dir,
  "Individual_PDF"
)
individual_svg_dir <- file.path(
  output_dir,
  "Individual_SVG"
)
combined_pdf_dir <- file.path(
  output_dir,
  "Combined_PDF"
)
combined_svg_dir <- file.path(
  output_dir,
  "Combined_SVG"
)
check_dir <- file.path(
  output_dir,
  "Matching_check"
)
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  individual_pdf_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  individual_svg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  combined_pdf_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  combined_svg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  check_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# ============================================================
# 3. 4
# ============================================================
focal_taxa <- c(
  "Staphylococcus capitis",
  "Pseudomonas syncyanea",
  "Luteitalea pratensis",
  "Paraflavitalea soli"
)
# ============================================================
# 4.
# ============================================================
outcome_variables <- c(
  "Yield",
  "RS_TN",
  "RS_TK",
  "RS_AP"
)
outcome_full_names <- c(
  "Yield" = "Yield",
  "RS_TN" =
    "Rhizosphere total nitrogen",
  "RS_TK" =
    "Rhizosphere total potassium",
  "RS_AP" =
    "Rhizosphere available phosphorus"
)
# ============================================================
# 5.
# ============================================================
color_control <- "#55B6A8"
color_mowing <- "#E47B6B"
color_dark <- "#222222"
color_ci <- "#BDBDBD"
# ============================================================
# 6.
# ============================================================
species_raw <- read_csv(
  species_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)
plant_raw <- read_excel(
  excel_file,
  sheet = "Plant",
  .name_repair = "minimal"
)
soil_raw <- read_excel(
  excel_file,
  sheet = "Soil",
  .name_repair = "minimal"
)
message(
  ":",
  nrow(species_raw),
  " × ",
  ncol(species_raw),
  ""
)
message(
  ":",
  nrow(plant_raw),
  " × ",
  ncol(plant_raw),
  ""
)
message(
  ":",
  nrow(soil_raw),
  " × ",
  ncol(soil_raw),
  ""
)
# ============================================================
# 7.
# ============================================================
required_species_columns <- c(
  "Species"
)
required_plant_columns <- c(
  "AnalysisID",
  "Tretment",
  "Year",
  "Yield_t_ha"
)
required_soil_columns <- c(
  "Soil Type",
  "Tretment",
  "Year",
  "Total Nitrogen g/kg",
  "Total potassium (g/kg)",
  "Available Phosphorus (mg/kg)"
)
missing_species_columns <- setdiff(
  required_species_columns,
  names(species_raw)
)
missing_plant_columns <- setdiff(
  required_plant_columns,
  names(plant_raw)
)
missing_soil_columns <- setdiff(
  required_soil_columns,
  names(soil_raw)
)
if (length(missing_species_columns) > 0) {
  stop(
    paste0(
      ":\n",
      paste(
        missing_species_columns,
        collapse = "\n"
      )
    )
  )
}
if (length(missing_plant_columns) > 0) {
  stop(
    paste0(
      ":\n",
      paste(
        missing_plant_columns,
        collapse = "\n"
      )
    )
  )
}
if (length(missing_soil_columns) > 0) {
  stop(
    paste0(
      ":\n",
      paste(
        missing_soil_columns,
        collapse = "\n"
      )
    )
  )
}
# ============================================================
# 8. 4
# ============================================================
missing_taxa <- setdiff(
  focal_taxa,
  as.character(
    species_raw$Species
  )
)
if (length(missing_taxa) > 0) {
  stop(
    paste0(
      ":\n",
      paste(
        missing_taxa,
        collapse = "\n"
      )
    )
  )
}
message(
  "4:",
  paste(
    focal_taxa,
    collapse = ";"
  )
)
# ============================================================
# 9. RT
# ============================================================
rt_sample_columns <- names(species_raw)[
  str_detect(
    names(species_raw),
    "^RT-[HG](CK|[0-9]+)-[0-9]+$"
  )
]
if (length(rt_sample_columns) == 0) {
  stop(
    paste0(
      "RT.\n",
      ":RT-H6-1,RT-G6-1,RT-HCK-1."
    )
  )
}
message(
  "RT:",
  length(rt_sample_columns)
)
# ============================================================
# 10.
# ============================================================
sample_total_check <- species_raw %>%
  summarise(
    across(
      all_of(
        rt_sample_columns
      ),
      ~sum(
        suppressWarnings(
          as.numeric(.x)
        ),
        na.rm = TRUE
      )
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "SampleID",
    values_to = "Total_abundance"
  )
write_csv(
  sample_total_check,
  file.path(
    check_dir,
    "RT_sample_abundance_total_check.csv"
  )
)
median_total_abundance <- median(
  sample_total_check$Total_abundance,
  na.rm = TRUE
)
message(
  "RT:",
  round(
    median_total_abundance,
    4
  )
)
if (
  median_total_abundance > 90 &&
  median_total_abundance < 110
) {
  abundance_multiplier <- 1
  message(
    ":,100."
  )
} else if (
  median_total_abundance > 0.9 &&
  median_total_abundance < 1.1
) {
  abundance_multiplier <- 100
  warning(
    paste0(
      ":0—1,",
      "100."
    )
  )
} else {
  abundance_multiplier <- 1
  warning(
    paste0(
      "100,1.\n",
      ",."
    )
  )
}
# ============================================================
# 11. 4
# ============================================================
sum_or_na <- function(x) {
  x <- suppressWarnings(
    as.numeric(x)
  )
  if (all(is.na(x))) {
    return(
      NA_real_
    )
  }
  sum(
    x,
    na.rm = TRUE
  )
}
species_four <- species_raw %>%
  filter(
    Species %in% focal_taxa
  ) %>%
  select(
    Species,
    all_of(
      rt_sample_columns
    )
  ) %>%
  group_by(
    Species
  ) %>%
  summarise(
    across(
      all_of(
        rt_sample_columns
      ),
      sum_or_na
    ),
    .groups = "drop"
  )
# ============================================================
# 12. RT
# ============================================================
rt_long <- species_four %>%
  pivot_longer(
    cols = all_of(
      rt_sample_columns
    ),
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  extract(
    col = SampleID,
    into = c(
      "Treatment",
      "Year_raw",
      "Rep"
    ),
    regex = "^RT-([HG])(?:CK|(\\d+))-(\\d+)$",
    remove = FALSE
  ) %>%
  mutate(
    Year = if_else(
      str_detect(
        SampleID,
        "CK"
      ),
      0L,
      suppressWarnings(
        as.integer(
          Year_raw
        )
      )
    ),
    Rep = suppressWarnings(
      as.integer(
        Rep
      )
    ),
    Treatment = as.character(
      Treatment
    ),
    Abundance = suppressWarnings(
      as.numeric(
        Abundance
      )
    ) * abundance_multiplier
  ) %>%
  select(
    SampleID,
    Treatment,
    Year,
    Rep,
    Species,
    Abundance
  )
# ============================================================
# 13. RT
# ============================================================
unparsed_rt_samples <- rt_long %>%
  filter(
    is.na(Treatment) |
      is.na(Year) |
      is.na(Rep)
  ) %>%
  distinct(
    SampleID
  )
write_csv(
  unparsed_rt_samples,
  file.path(
    check_dir,
    "Unparsed_RT_sample_IDs.csv"
  )
)
if (nrow(unparsed_rt_samples) > 0) {
  stop(
    paste0(
      "RT.\n",
      ":Unparsed_RT_sample_IDs.csv"
    )
  )
}
# ============================================================
# 14. RT
# ============================================================
rt_wide <- rt_long %>%
  pivot_wider(
    names_from = Species,
    values_from = Abundance
  )
# RT
rt_duplicate_keys <- rt_wide %>%
  count(
    Treatment,
    Year,
    Rep,
    name = "Duplicate_number"
  ) %>%
  filter(
    Duplicate_number > 1
  )
write_csv(
  rt_duplicate_keys,
  file.path(
    check_dir,
    "Duplicate_RT_matching_keys.csv"
  )
)
if (nrow(rt_duplicate_keys) > 0) {
  stop(
    paste0(
      "RTTreatment-Year-Rep.\n",
      ":Duplicate_RT_matching_keys.csv"
    )
  )
}
message(
  "RT:",
  nrow(rt_wide)
)
# ============================================================
# 15.
# ============================================================
plant_data <- plant_raw %>%
  transmute(
    AnalysisID = as.character(
      ``
    ),
    Treatment = str_sub(
      str_trim(
        as.character(
          Tretment
        )
      ),
      1,
      1
    ),
    Year = suppressWarnings(
      as.integer(
        Year
      )
    ),
    Rep = suppressWarnings(
      as.integer(
        str_extract(
          as.character(
            ``
          ),
          "\\d+$"
        )
      )
    ),
    Yield = suppressWarnings(
      as.numeric(
        `Yield_t_ha`
      )
    )
  ) %>%
  filter(
    Treatment %in% c(
      "G",
      "H"
    ),
    !is.na(Year),
    !is.na(Rep),
    is.finite(
      Yield
    )
  )
plant_duplicate_keys <- plant_data %>%
  count(
    Treatment,
    Year,
    Rep,
    name = "Duplicate_number"
  ) %>%
  filter(
    Duplicate_number > 1
  )
write_csv(
  plant_duplicate_keys,
  file.path(
    check_dir,
    "Duplicate_plant_matching_keys.csv"
  )
)
if (nrow(plant_duplicate_keys) > 0) {
  stop(
    paste0(
      "Treatment-Year-Rep.\n",
      ":Duplicate_plant_matching_keys.csv"
    )
  )
}
message(
  ":",
  nrow(plant_data)
)
# ============================================================
# 16. RS
# ============================================================
# ExcelRep.
# Soil Type × Treatment × Year,
# Rep 1—5.
# Excel51—5.
soil_group_check <- soil_raw %>%
  mutate(
    Treatment = str_sub(
      str_trim(
        as.character(
          Tretment
        )
      ),
      1,
      1
    ),
    Year = suppressWarnings(
      as.integer(
        Year
      )
    )
  ) %>%
  count(
    `Soil Type`,
    Treatment,
    Year,
    name = "Replicate_number"
  )
write_csv(
  soil_group_check,
  file.path(
    check_dir,
    "Soil_group_replicate_check.csv"
  )
)
invalid_soil_groups <- soil_group_check %>%
  filter(
    Replicate_number != 5
  )
if (nrow(invalid_soil_groups) > 0) {
  warning(
    paste0(
      "5.\n",
      ":Soil_group_replicate_check.csv"
    )
  )
}
rs_nutrients <- soil_raw %>%
  mutate(
    Original_row = row_number(),
    Treatment = str_sub(
      str_trim(
        as.character(
          Tretment
        )
      ),
      1,
      1
    ),
    Year = suppressWarnings(
      as.integer(
        Year
      )
    )
  ) %>%
  filter(
    `Soil Type` == "RS",
    Treatment %in% c(
      "G",
      "H"
    ),
    !is.na(Year)
  ) %>%
  group_by(
    `Soil Type`,
    Treatment,
    Year
  ) %>%
  arrange(
    Original_row,
    .by_group = TRUE
  ) %>%
  mutate(
    Rep = row_number()
  ) %>%
  ungroup() %>%
  transmute(
    Treatment,
    Year,
    Rep,
    RS_TN = suppressWarnings(
      as.numeric(
        `Total Nitrogen g/kg`
      )
    ),
    RS_TK = suppressWarnings(
      as.numeric(
        `Total potassium (g/kg)`
      )
    ),
    RS_AP = suppressWarnings(
      as.numeric(
        `Available Phosphorus (mg/kg)`
      )
    )
  )
# RS
soil_duplicate_keys <- rs_nutrients %>%
  count(
    Treatment,
    Year,
    Rep,
    name = "Duplicate_number"
  ) %>%
  filter(
    Duplicate_number > 1
  )
write_csv(
  soil_duplicate_keys,
  file.path(
    check_dir,
    "Duplicate_RS_nutrient_matching_keys.csv"
  )
)
if (nrow(soil_duplicate_keys) > 0) {
  stop(
    paste0(
      "RS.\n",
      ":Duplicate_RS_nutrient_matching_keys.csv"
    )
  )
}
message(
  "RS:",
  nrow(rs_nutrients)
)
# ============================================================
# 17. RT,RS
# ============================================================
matching_keys <- c(
  "Treatment",
  "Year",
  "Rep"
)
matched_data <- rt_wide %>%
  inner_join(
    plant_data,
    by = matching_keys
  ) %>%
  inner_join(
    rs_nutrients,
    by = matching_keys
  )
if (nrow(matched_data) == 0) {
  stop(
    paste0(
      "RT,RS.\n",
      "Treatment,YearRep."
    )
  )
}
# ============================================================
# 18.
# ============================================================
unmatched_rt_vs_plant <- rt_wide %>%
  anti_join(
    plant_data,
    by = matching_keys
  )
unmatched_rt_vs_soil <- rt_wide %>%
  anti_join(
    rs_nutrients,
    by = matching_keys
  )
unmatched_plant_vs_rt <- plant_data %>%
  anti_join(
    rt_wide,
    by = matching_keys
  )
unmatched_soil_vs_rt <- rs_nutrients %>%
  anti_join(
    rt_wide,
    by = matching_keys
  )
write_csv(
  unmatched_rt_vs_plant,
  file.path(
    check_dir,
    "Unmatched_RT_vs_plant.csv"
  )
)
write_csv(
  unmatched_rt_vs_soil,
  file.path(
    check_dir,
    "Unmatched_RT_vs_RS_nutrients.csv"
  )
)
write_csv(
  unmatched_plant_vs_rt,
  file.path(
    check_dir,
    "Unmatched_plant_vs_RT.csv"
  )
)
write_csv(
  unmatched_soil_vs_rt,
  file.path(
    check_dir,
    "Unmatched_RS_nutrients_vs_RT.csv"
  )
)
write_csv(
  matched_data,
  file.path(
    output_dir,
    "Four_RT_taxa_Yield_RS_nutrients_matched.csv"
  )
)
message(
  "RT:",
  nrow(rt_wide)
)
message(
  ":",
  nrow(matched_data)
)
message(
  "RT:",
  nrow(unmatched_rt_vs_plant)
)
message(
  "RTRS:",
  nrow(unmatched_rt_vs_soil)
)
# ============================================================
# 19.
# ============================================================
matched_data <- matched_data %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c(
        "G",
        "H"
      )
    ),
    Year = suppressWarnings(
      as.integer(
        Year
      )
    ),
    Rep = suppressWarnings(
      as.integer(
        Rep
      )
    ),
    across(
      all_of(
        c(
          focal_taxa,
          outcome_variables
        )
      ),
      ~suppressWarnings(
        as.numeric(
          .x
        )
      )
    )
  )
# ============================================================
# 20.
# ============================================================
variable_check <- tibble(
  Variable = c(
    focal_taxa,
    outcome_variables
  ),
  Valid_n = map_int(
    c(
      focal_taxa,
      outcome_variables
    ),
    function(variable_name) {
      sum(
        is.finite(
          matched_data[[variable_name]]
        )
      )
    }
  ),
  Unique_n = map_int(
    c(
      focal_taxa,
      outcome_variables
    ),
    function(variable_name) {
      valid_values <- matched_data[[variable_name]][
        is.finite(
          matched_data[[variable_name]]
        )
      ]
      n_distinct(
        valid_values
      )
    }
  )
)
write_csv(
  variable_check,
  file.path(
    output_dir,
    "Variable_data_check.csv"
  )
)
# ============================================================
# 21.
# ============================================================
run_linear_regression <- function(
    data,
    taxon_name,
    outcome_name
) {
  regression_data <- data %>%
    transmute(
      Abundance = suppressWarnings(
        as.numeric(
          .data[[taxon_name]]
        )
      ),
      Response = suppressWarnings(
        as.numeric(
          .data[[outcome_name]]
        )
      )
    ) %>%
    filter(
      is.finite(
        Abundance
      ),
      is.finite(
        Response
      )
    )
  sample_number <- nrow(
    regression_data
  )
  if (
    sample_number < 3 ||
    n_distinct(
      regression_data$Abundance
    ) < 2 ||
    n_distinct(
      regression_data$Response
    ) < 2
  ) {
    return(
      tibble(
        Taxon = taxon_name,
        Outcome = outcome_name,
        n = sample_number,
        Intercept = NA_real_,
        Slope = NA_real_,
        Slope_SE = NA_real_,
        Slope_lower_95CI = NA_real_,
        Slope_upper_95CI = NA_real_,
        R2 = NA_real_,
        Adjusted_R2 = NA_real_,
        P_value = NA_real_
      )
    )
  }
  regression_model <- lm(
    Response ~ Abundance,
    data = regression_data
  )
  regression_summary <- summary(
    regression_model
  )
  slope_confidence_interval <- confint(
    regression_model,
    parm = "Abundance",
    level = 0.95
  )
  tibble(
    Taxon = taxon_name,
    Outcome = outcome_name,
    n = sample_number,
    Intercept = unname(
      coef(
        regression_model
      )[["(Intercept)"]]
    ),
    Slope = unname(
      coef(
        regression_model
      )[["Abundance"]]
    ),
    Slope_SE = regression_summary$coefficients[
      "Abundance",
      "Std. Error"
    ],
    Slope_lower_95CI =
      unname(
        slope_confidence_interval[
          1,
          1
        ]
      ),
    Slope_upper_95CI =
      unname(
        slope_confidence_interval[
          1,
          2
        ]
      ),
    R2 = regression_summary$r.squared,
    Adjusted_R2 =
      regression_summary$adj.r.squared,
    P_value =
      regression_summary$coefficients[
        "Abundance",
        "Pr(>|t|)"
      ]
  )
}
# ============================================================
# 22. 16
# ============================================================
regression_combinations <- expand_grid(
  Outcome = outcome_variables,
  Taxon = focal_taxa
)
regression_results <- map2_dfr(
  regression_combinations$Taxon,
  regression_combinations$Outcome,
  function(
    taxon_name,
    outcome_name
  ) {
    run_linear_regression(
      data = matched_data,
      taxon_name = taxon_name,
      outcome_name = outcome_name
    )
  }
)
# ============================================================
# 23. FDR
# ============================================================
# FDR_all_16
# 16BH.
# FDR_within_outcome
# ,4BH.
regression_results <- regression_results %>%
  mutate(
    FDR_all_16 = p.adjust(
      P_value,
      method = "BH"
    )
  ) %>%
  group_by(
    Outcome
  ) %>%
  mutate(
    FDR_within_outcome = p.adjust(
      P_value,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    Direction = case_when(
      is.na(Slope) ~ NA_character_,
      Slope > 0 ~ "Positive",
      Slope < 0 ~ "Negative",
      TRUE ~ "No direction"
    ),
    Significant_all_16 = case_when(
      is.na(FDR_all_16) ~ "NA",
      FDR_all_16 < 0.001 ~ "***",
      FDR_all_16 < 0.01 ~ "**",
      FDR_all_16 < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    Outcome_full_name =
      unname(
        outcome_full_names[
          Outcome
        ]
      )
  ) %>%
  arrange(
    factor(
      Outcome,
      levels = outcome_variables
    ),
    FDR_all_16,
    desc(
      abs(
        Slope
      )
    )
  )
write_csv(
  regression_results,
  file.path(
    output_dir,
    "Four_taxa_16_linear_regression_results.csv"
  )
)
significant_results <- regression_results %>%
  filter(
    !is.na(
      FDR_all_16
    ),
    FDR_all_16 < 0.05
  )
write_csv(
  significant_results,
  file.path(
    output_dir,
    "Four_taxa_significant_regressions_FDR_all_16_lt_0.05.csv"
  )
)
message(
  ":",
  nrow(regression_results)
)
message(
  "16,FDR < 0.05:",
  nrow(significant_results)
)
# ============================================================
# 24.
# ============================================================
format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return(
      "P = NA"
    )
  }
  if (p_value < 0.001) {
    return(
      "P < 0.001"
    )
  }
  paste0(
    "P = ",
    sprintf(
      "%.3f",
      p_value
    )
  )
}
format_fdr_value <- function(fdr_value) {
  if (is.na(fdr_value)) {
    return(
      "FDR = NA"
    )
  }
  if (fdr_value < 0.001) {
    return(
      "FDR < 0.001"
    )
  }
  paste0(
    "FDR = ",
    sprintf(
      "%.3f",
      fdr_value
    )
  )
}
safe_filename <- function(x) {
  x %>%
    str_replace_all(
      "\\s+",
      "_"
    ) %>%
    str_replace_all(
      "[^A-Za-z0-9_.-]",
      "_"
    ) %>%
    str_replace_all(
      "_+",
      "_"
    )
}
# ============================================================
# 25. Y
# ============================================================
get_y_axis_label <- function(outcome_name) {
  switch(
    outcome_name,
    "Yield" = expression(
      Yield~(t~ha^{-1})
    ),
    "RS_TN" = expression(
      "Rhizosphere TN"~(g~kg^{-1})
    ),
    "RS_TK" = expression(
      "Rhizosphere TK"~(g~kg^{-1})
    ),
    "RS_AP" = expression(
      "Rhizosphere AP"~(mg~kg^{-1})
    ),
    outcome_name
  )
}
# ============================================================
# 26.
# ============================================================
theme_regression <- theme_classic(
  base_family = "serif"
) +
  theme(
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    panel.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    panel.border = element_rect(
      fill = NA,
      color = color_dark,
      linewidth = 0.55
    ),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(
      color = color_dark,
      size = 8
    ),
    axis.title = element_text(
      color = color_dark,
      size = 9
    ),
    axis.title.x = element_text(
      margin = margin(
        t = 6
      )
    ),
    axis.title.y = element_text(
      margin = margin(
        r = 6
      )
    ),
    plot.title = element_text(
      color = color_dark,
      size = 10,
      hjust = 0.5,
      face = "italic",
      margin = margin(
        b = 5
      )
    ),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(
      color = color_dark,
      size = 8
    ),
    legend.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    legend.box.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    legend.key = element_rect(
      fill = "transparent",
      color = NA
    ),
    plot.tag = element_text(
      color = color_dark,
      size = 11,
      face = "bold"
    ),
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )
# ============================================================
# 27.
# ============================================================
make_regression_plot <- function(
    data,
    taxon_name,
    outcome_name
) {
  plot_data <- data %>%
    transmute(
      Treatment = factor(
        as.character(
          Treatment
        ),
        levels = c(
          "G",
          "H"
        )
      ),
      Abundance = suppressWarnings(
        as.numeric(
          .data[[taxon_name]]
        )
      ),
      Response = suppressWarnings(
        as.numeric(
          .data[[outcome_name]]
        )
      )
    ) %>%
    filter(
      is.finite(
        Abundance
      ),
      is.finite(
        Response
      ),
      !is.na(
        Treatment
      )
    )
  stat_row <- regression_results %>%
    filter(
      Taxon == taxon_name,
      Outcome == outcome_name
    ) %>%
    slice(
      1
    )
  stat_text <- paste0(
    "R² = ",
    ifelse(
      is.na(
        stat_row$R2[[1]]
      ),
      "NA",
      sprintf(
        "%.3f",
        stat_row$R2[[1]]
      )
    ),
    "\n",
    format_p_value(
      stat_row$P_value[[1]]
    ),
    "\n",
    format_fdr_value(
      stat_row$FDR_all_16[[1]]
    )
  )
  ggplot(
    plot_data,
    aes(
      x = Abundance,
      y = Response
    )
  ) +
 # 95%
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      level = 0.95,
      color = color_dark,
      fill = color_ci,
      alpha = 0.32,
      linewidth = 0.70
    ) +
 # GH
    geom_point(
      aes(
        color = Treatment
      ),
      size = 2.3,
      alpha = 0.86,
      shape = 16
    ) +
    scale_color_manual(
      values = c(
        "G" = color_control,
        "H" = color_mowing
      ),
      breaks = c(
        "G",
        "H"
      ),
      labels = c(
        "Control",
        "Mowing"
      ),
      drop = FALSE
    ) +
    annotate(
      geom = "text",
      x = Inf,
      y = Inf,
      label = stat_text,
      hjust = 1.10,
      vjust = 1.12,
      color = color_dark,
      size = 3.0,
      lineheight = 0.95,
      family = "serif"
    ) +
    labs(
      title = taxon_name,
      x = "Relative abundance (%)",
      y = get_y_axis_label(
        outcome_name
      )
    ) +
    scale_x_continuous(
      expand = expansion(
        mult = c(
          0.05,
          0.18
        )
      )
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.06,
          0.22
        )
      )
    ) +
    theme_regression
}
# ============================================================
# 28. 164
# ============================================================
for (outcome_name in outcome_variables) {
  outcome_pdf_dir <- file.path(
    individual_pdf_dir,
    outcome_name
  )
  outcome_svg_dir <- file.path(
    individual_svg_dir,
    outcome_name
  )
  dir.create(
    outcome_pdf_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  dir.create(
    outcome_svg_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  one_outcome_plots <- list()
  for (taxon_name in focal_taxa) {
    p_single <- make_regression_plot(
      data = matched_data,
      taxon_name = taxon_name,
      outcome_name = outcome_name
    )
    one_outcome_plots[[
      taxon_name
    ]] <- p_single
    file_stub <- safe_filename(
      taxon_name
    )
 # --------------------------------------------------------
 # PDF
 # --------------------------------------------------------
    pdf_file <- file.path(
      outcome_pdf_dir,
      paste0(
        file_stub,
        "_vs_",
        outcome_name,
        "_linear_regression.pdf"
      )
    )
    ggsave(
      filename = pdf_file,
      plot = p_single,
      width = 4.2,
      height = 3.8,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = "transparent"
    )
 # --------------------------------------------------------
 # SVG
 # --------------------------------------------------------
    svg_file <- file.path(
      outcome_svg_dir,
      paste0(
        file_stub,
        "_vs_",
        outcome_name,
        "_linear_regression.svg"
      )
    )
    grDevices::svg(
      filename = svg_file,
      width = 4.2,
      height = 3.8,
      bg = "transparent",
      onefile = TRUE
    )
    print(
      p_single
    )
    grDevices::dev.off()
    message(
      ":",
      taxon_name,
      " × ",
      outcome_name
    )
  }
 # ==========================================================
 # 42×2
 # ==========================================================
  combined_plot <- wrap_plots(
    one_outcome_plots,
    ncol = 2,
    guides = "collect"
  ) +
    plot_annotation(
      tag_levels = "A"
    ) &
    theme(
      legend.position = "bottom"
    )
 # PDF
  combined_pdf_file <- file.path(
    combined_pdf_dir,
    paste0(
      "Four_taxa_vs_",
      outcome_name,
      "_combined.pdf"
    )
  )
  ggsave(
    filename = combined_pdf_file,
    plot = combined_plot,
    width = 8.5,
    height = 7.2,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "transparent"
  )
 # SVG
  combined_svg_file <- file.path(
    combined_svg_dir,
    paste0(
      "Four_taxa_vs_",
      outcome_name,
      "_combined.svg"
    )
  )
  grDevices::svg(
    filename = combined_svg_file,
    width = 8.5,
    height = 7.2,
    bg = "transparent",
    onefile = TRUE
  )
  print(
    combined_plot
  )
  grDevices::dev.off()
}
# ============================================================
# 29.
# ============================================================
regression_summary_wide <- regression_results %>%
  transmute(
    Taxon,
    Outcome,
    Result = paste0(
      "Slope = ",
      ifelse(
        is.na(Slope),
        "NA",
        sprintf(
          "%.5f",
          Slope
        )
      ),
      "; R2 = ",
      ifelse(
        is.na(R2),
        "NA",
        sprintf(
          "%.3f",
          R2
        )
      ),
      "; P = ",
      ifelse(
        is.na(P_value),
        "NA",
        format(
          P_value,
          digits = 3,
          scientific = TRUE
        )
      ),
      "; FDR_all_16 = ",
      ifelse(
        is.na(FDR_all_16),
        "NA",
        format(
          FDR_all_16,
          digits = 3,
          scientific = TRUE
        )
      ),
      "; FDR_within_outcome = ",
      ifelse(
        is.na(FDR_within_outcome),
        "NA",
        format(
          FDR_within_outcome,
          digits = 3,
          scientific = TRUE
        )
      )
    )
  ) %>%
  pivot_wider(
    names_from = Outcome,
    values_from = Result
  )
write_csv(
  regression_summary_wide,
  file.path(
    output_dir,
    "Four_taxa_regression_summary_wide.csv"
  )
)
# ============================================================
# 30.
# ============================================================
message(
  "============================================================"
)
message(
  "."
)
message(
  ":",
  nrow(matched_data)
)
message(
  ":",
  nrow(regression_results)
)
message(
  "PDF:",
  normalizePath(
    individual_pdf_dir
  )
)
message(
  "SVG:",
  normalizePath(
    individual_svg_dir
  )
)
message(
  "PDF:",
  normalizePath(
    combined_pdf_dir
  )
)
message(
  "SVG:",
  normalizePath(
    combined_svg_dir
  )
)
message(
  ":",
  normalizePath(
    file.path(
      output_dir,
      "Four_taxa_16_linear_regression_results.csv"
    )
  )
)
message(
  ":",
  normalizePath(
    file.path(
      output_dir,
      "Four_RT_taxa_Yield_RS_nutrients_matched.csv"
    )
  )
)
message(
  "============================================================"
)
