# -*- coding: UTF-8 -*-
# ============================================================
# = RS
# = 4RT
# = BS
# Three-column association network
# Left = high-importance RS indicators
# Middle = four yield-associated RT taxa
# Right = high-importance BS indicators
# 4RS/BSSpearman
# 40BHFDR < 0.05.
# Edge rule
# Spearman associations between four taxa and RS/BS indicators
# only associations with global BH-adjusted FDR < 0.05 are shown.
# = ; =
# = |Spearman rho|
# Red edges = positive associations
# Blue edges = negative associations
# Edge width = |Spearman rho|
# ============================================================
# ============================================================
# 0. R / Required packages
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
  "scales"
)
missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "R:\n",
      "install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
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
  library(scales)
})
# ============================================================
# 1. / File paths
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
output_dir <- file.path(
  work_dir,
  "RS_taxa_BS_three_column_network"
)
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!file.exists(excel_file)) {
  stop("Input file not found: main_plant_soil_data.xlsx")
}
if (!file.exists(species_file)) {
  stop("Input file not found: species_abundance.csv")
}
# ============================================================
# 2. 4
# Four yield-associated taxa
# ============================================================
focal_taxa <- c(
  "Staphylococcus capitis",
  "Pseudomonas syncyanea",
  "Luteitalea pratensis",
  "Paraflavitalea soli"
)
# ============================================================
# 3.
# Indicators and contribution-based order
# ============================================================
# RS
# RS: ordered from higher to lower contribution
rs_indicator_order <- c(
  "TK",
  "Turnover",
  "AP",
  "TN",
  "MBC:MBN"
)
# BS:,
# FDR.
# BS: high-importance indicators retained when at least one
# FDR-significant association with the four taxa is present.
# AK,4FDR < 0.05,
# ,.
bs_indicator_order <- c(
  "pH",
  "Turnover",
  "AP",
  "TK",
  "MBC:MBN"
)
# AK,
# bs_indicator_order <- c("pH", "Turnover", "AP", "AK", "TK")
# ============================================================
# 4. / Statistical settings
# ============================================================
fdr_threshold <- 0.05
# TRUE:FDR
# FALSE
show_only_significant_edges <- TRUE
# ============================================================
# 5. / Read data
# ============================================================
species_raw <- read_csv(
  species_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)
soil_raw <- read_excel(
  excel_file,
  sheet = "Soil",
  .name_repair = "minimal"
)
carbon_raw <- read_excel(
  excel_file,
  sheet = "CarbonUse",
  .name_repair = "minimal"
)
# ============================================================
# 6. / Check required data
# ============================================================
if (!"Species" %in% names(species_raw)) {
  stop("SpeciesSpecies.")
}
missing_taxa <- setdiff(
  focal_taxa,
  as.character(species_raw$Species)
)
if (length(missing_taxa) > 0) {
  stop(
    paste0(
      "Species:",
      paste(missing_taxa, collapse = ";")
    )
  )
}
required_soil_columns <- c(
  "Soil Type",
  "Tretment",
  "Year",
  "Total Nitrogen g/kg",
  "Total potassium (g/kg)",
  "Available Phosphorus (mg/kg)",
  "Available Potassium (mg/kg)",
  "PH"
)
required_carbon_columns <- c(
  "Soil Type",
  "Tretment",
  "Year",
  "Microbial biomass C/N",
  "Microbial respiration rate (ng C g-1 soil h-1)",
  "Turnover rate"
)
missing_soil_columns <- setdiff(
  required_soil_columns,
  names(soil_raw)
)
missing_carbon_columns <- setdiff(
  required_carbon_columns,
  names(carbon_raw)
)
if (length(missing_soil_columns) > 0) {
  stop(
    paste0(
      ":",
      paste(missing_soil_columns, collapse = ";")
    )
  )
}
if (length(missing_carbon_columns) > 0) {
  stop(
    paste0(
      ":",
      paste(missing_carbon_columns, collapse = ";")
    )
  )
}
# ============================================================
# 7. 4RT
# Extract RT abundances of the four taxa
# ============================================================
rt_sample_columns <- names(species_raw)[
  str_detect(
    names(species_raw),
    "^RT-[HG](CK|[0-9]+)-[0-9]+$"
  )
]
if (length(rt_sample_columns) == 0) {
  stop("RT.")
}
sum_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) {
    return(NA_real_)
  }
  sum(x, na.rm = TRUE)
}
species_four <- species_raw %>%
  filter(
    Species %in% focal_taxa
  ) %>%
  select(
    Species,
    all_of(rt_sample_columns)
  ) %>%
  group_by(
    Species
  ) %>%
  summarise(
    across(
      all_of(rt_sample_columns),
      sum_or_na
    ),
    .groups = "drop"
  )
rt_long <- species_four %>%
  pivot_longer(
    cols = all_of(rt_sample_columns),
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
    Treatment = as.character(Treatment),
    Year = if_else(
      str_detect(SampleID, "CK"),
      0L,
      suppressWarnings(as.integer(Year_raw))
    ),
    Rep = suppressWarnings(
      as.integer(Rep)
    ),
    Abundance = suppressWarnings(
      as.numeric(Abundance)
    )
  ) %>%
  select(
    SampleID,
    Treatment,
    Year,
    Rep,
    Species,
    Abundance
  )
if (
  any(
    is.na(rt_long$Treatment) |
      is.na(rt_long$Year) |
      is.na(rt_long$Rep)
  )
) {
  stop("RT,.")
}
rt_wide <- rt_long %>%
  pivot_wider(
    names_from = Species,
    values_from = Abundance
  )
# ============================================================
# 8.
# Prepare soil and microbial indicators
# ============================================================
prepare_compartment_indicators <- function(
    compartment_name
) {
  soil_part <- soil_raw %>%
    mutate(
      Original_row = row_number(),
      Treatment = str_sub(
        str_trim(as.character(Tretment)),
        1,
        1
      ),
      Year = suppressWarnings(
        as.integer(Year)
      )
    ) %>%
    filter(
      `Soil Type` == compartment_name,
      Treatment %in% c("G", "H"),
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
      TN = suppressWarnings(
        as.numeric(`Total Nitrogen g/kg`)
      ),
      TK = suppressWarnings(
        as.numeric(`Total potassium (g/kg)`)
      ),
      AP = suppressWarnings(
        as.numeric(`Available Phosphorus (mg/kg)`)
      ),
      AK = suppressWarnings(
        as.numeric(`Available Potassium (mg/kg)`)
      ),
      pH = suppressWarnings(
        as.numeric(PH)
      )
    )
  carbon_part <- carbon_raw %>%
    mutate(
      Original_row = row_number(),
      Treatment = str_sub(
        str_trim(as.character(Tretment)),
        1,
        1
      ),
      Year = suppressWarnings(
        as.integer(Year)
      )
    ) %>%
    filter(
      `Soil Type` == compartment_name,
      Treatment %in% c("G", "H"),
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
      `MBC:MBN` = suppressWarnings(
        as.numeric(`Microbial biomass C/N`)
      ),
      BR = suppressWarnings(
        as.numeric(
          `Microbial respiration rate (ng C g-1 soil h-1)`
        )
      ),
      Turnover = suppressWarnings(
        as.numeric(`Turnover rate`)
      )
    )
  soil_part %>%
    inner_join(
      carbon_part,
      by = c(
        "Treatment",
        "Year",
        "Rep"
      )
    )
}
rs_indicators <- prepare_compartment_indicators(
  "RS"
)
bs_indicators <- prepare_compartment_indicators(
  "BS"
)
# ============================================================
# 9. RTRS/BS
# Match RT taxa with RS and BS indicators
# ============================================================
matching_keys <- c(
  "Treatment",
  "Year",
  "Rep"
)
rt_rs_data <- rt_wide %>%
  inner_join(
    rs_indicators,
    by = matching_keys
  )
rt_bs_data <- rt_wide %>%
  inner_join(
    bs_indicators,
    by = matching_keys
  )
if (nrow(rt_rs_data) == 0) {
  stop("RTRS.")
}
if (nrow(rt_bs_data) == 0) {
  stop("RTBS.")
}
message(
  "RT–RS:",
  nrow(rt_rs_data)
)
message(
  "RT–BS:",
  nrow(rt_bs_data)
)
# ============================================================
# 10. Spearman
# Spearman correlation function
# ============================================================
calculate_one_correlation <- function(
    data,
    compartment_name,
    taxon_name,
    indicator_name
) {
  x <- suppressWarnings(
    as.numeric(data[[taxon_name]])
  )
  y <- suppressWarnings(
    as.numeric(data[[indicator_name]])
  )
  valid <- complete.cases(x, y) &
    is.finite(x) &
    is.finite(y)
  x <- x[valid]
  y <- y[valid]
  if (
    length(x) < 3 ||
      n_distinct(x) < 2 ||
      n_distinct(y) < 2
  ) {
    return(
      tibble(
        Compartment = compartment_name,
        Taxon = taxon_name,
        Indicator = indicator_name,
        rho = NA_real_,
        p = NA_real_,
        n = length(x)
      )
    )
  }
  test_result <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  tibble(
    Compartment = compartment_name,
    Taxon = taxon_name,
    Indicator = indicator_name,
    rho = unname(
      test_result$estimate
    ),
    p = test_result$p.value,
    n = length(x)
  )
}
# ============================================================
# 11. 40
# Calculate all 40 associations
# ============================================================
rs_combinations <- crossing(
  Taxon = focal_taxa,
  Indicator = rs_indicator_order
)
bs_combinations <- crossing(
  Taxon = focal_taxa,
  Indicator = bs_indicator_order
)
rs_results <- map2_dfr(
  rs_combinations$Taxon,
  rs_combinations$Indicator,
  function(taxon_name, indicator_name) {
    calculate_one_correlation(
      data = rt_rs_data,
      compartment_name = "RS",
      taxon_name = taxon_name,
      indicator_name = indicator_name
    )
  }
)
bs_results <- map2_dfr(
  bs_combinations$Taxon,
  bs_combinations$Indicator,
  function(taxon_name, indicator_name) {
    calculate_one_correlation(
      data = rt_bs_data,
      compartment_name = "BS",
      taxon_name = taxon_name,
      indicator_name = indicator_name
    )
  }
)
correlation_results <- bind_rows(
  rs_results,
  bs_results
) %>%
  mutate(
 # 40BH
    FDR_global = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  group_by(
    Compartment
  ) %>%
  mutate(
 # ,
    FDR_within_compartment = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    Significant = !is.na(FDR_global) &
      FDR_global < fdr_threshold,
    Direction = case_when(
      is.na(rho) ~ NA_character_,
      rho > 0 ~ "Positive",
      rho < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    Abs_rho = abs(rho)
  ) %>%
  arrange(
    Compartment,
    FDR_global,
    desc(Abs_rho)
  )
write_csv(
  correlation_results,
  file.path(
    output_dir,
    "RS_BS_four_taxa_Spearman_all_results.csv"
  )
)
significant_results <- correlation_results %>%
  filter(
    Significant
  )
write_csv(
  significant_results,
  file.path(
    output_dir,
    "RS_BS_four_taxa_Spearman_FDR_lt_0.05.csv"
  )
)
message(
  ":",
  nrow(correlation_results)
)
message(
  "FDR < 0.05:",
  nrow(significant_results)
)
# ============================================================
# 12.
# Build fixed three-column nodes
# ============================================================
# .
# Nodes are ordered from top to bottom by contribution rank.
rs_nodes <- tibble(
  Node = paste0(
    "RS__",
    rs_indicator_order
  ),
  Label = rs_indicator_order,
  Group = "RS",
  x = 1,
  y = rev(
    seq_along(
      rs_indicator_order
    )
  ),
  Rank = seq_along(
    rs_indicator_order
  )
)
taxa_nodes <- tibble(
  Node = paste0(
    "Taxon__",
    focal_taxa
  ),
  Label = focal_taxa,
  Group = "Taxon",
  x = 2.5,
  y = seq(
    from = 4.4,
    to = 1.6,
    length.out = length(focal_taxa)
  ),
  Rank = NA_integer_
)
bs_nodes <- tibble(
  Node = paste0(
    "BS__",
    bs_indicator_order
  ),
  Label = bs_indicator_order,
  Group = "BS",
  x = 4,
  y = rev(
    seq_along(
      bs_indicator_order
    )
  ),
  Rank = seq_along(
    bs_indicator_order
  )
)
nodes <- bind_rows(
  rs_nodes,
  taxa_nodes,
  bs_nodes
)
# ============================================================
# 13.
# Build edge table
# ============================================================
edge_source <- if (
  show_only_significant_edges
) {
  significant_results
} else {
  correlation_results %>%
    filter(
      !is.na(rho)
    )
}
edges <- edge_source %>%
  mutate(
    From = case_when(
      Compartment == "RS" ~
        paste0("RS__", Indicator),
      Compartment == "BS" ~
        paste0("Taxon__", Taxon)
    ),
    To = case_when(
      Compartment == "RS" ~
        paste0("Taxon__", Taxon),
      Compartment == "BS" ~
        paste0("BS__", Indicator)
    )
  ) %>%
  left_join(
    nodes %>%
      select(
        From = Node,
        x_from = x,
        y_from = y
      ),
    by = "From"
  ) %>%
  left_join(
    nodes %>%
      select(
        To = Node,
        x_to = x,
        y_to = y
      ),
    by = "To"
  ) %>%
  mutate(
    Edge_direction = factor(
      Direction,
      levels = c(
        "Positive",
        "Negative"
      )
    )
  )
write_csv(
  edges,
  file.path(
    output_dir,
    "Network_edges_used_for_plot.csv"
  )
)
# ============================================================
# 14.
# Plot colours and style
# ============================================================
color_rs <- "#C77C2E"
color_taxa <- "#2F7F9D"
color_bs <- "#7F7F7F"
color_positive <- "#C44E52"
color_negative <- "#3E6FB0"
color_border <- "#222222"
color_panel_fill <- "#FFFFFF"
# ============================================================
# 15.
# Draw three-column network
# ============================================================
network_plot <- ggplot() +
 # / Three column boxes
  annotate(
    geom = "rect",
    xmin = 0.45,
    xmax = 1.55,
    ymin = 0.45,
    ymax = 5.75,
    fill = color_panel_fill,
    color = color_border,
    linewidth = 0.8
  ) +
  annotate(
    geom = "rect",
    xmin = 1.80,
    xmax = 3.20,
    ymin = 0.45,
    ymax = 5.75,
    fill = color_panel_fill,
    color = color_border,
    linewidth = 0.8
  ) +
  annotate(
    geom = "rect",
    xmin = 3.45,
    xmax = 4.55,
    ymin = 0.45,
    ymax = 5.75,
    fill = color_panel_fill,
    color = color_border,
    linewidth = 0.8
  ) +
 # ,
 # Draw edges before nodes
  geom_segment(
    data = edges,
    aes(
      x = x_from,
      y = y_from,
      xend = x_to,
      yend = y_to,
      color = Edge_direction,
      linewidth = Abs_rho
    ),
    alpha = 0.56,
    lineend = "round"
  ) +
  scale_color_manual(
    name = "Association",
    values = c(
      "Positive" = color_positive,
      "Negative" = color_negative
    ),
    labels = c(
      "Positive",
      "Negative"
    ),
    drop = FALSE
  ) +
  scale_linewidth_continuous(
    name = expression("|Spearman "*rho*"|"),
    range = c(
      0.35,
      2.1
    ),
    limits = c(
      0,
      1
    )
  ) +
 # RS / RS nodes
  geom_point(
    data = rs_nodes,
    aes(
      x = x,
      y = y
    ),
    shape = 21,
    size = 6.0,
    stroke = 0.45,
    fill = color_rs,
    color = color_border
  ) +
 # / Middle taxon nodes
  geom_point(
    data = taxa_nodes,
    aes(
      x = x,
      y = y
    ),
    shape = 21,
    size = 7.0,
    stroke = 0.50,
    fill = color_taxa,
    color = color_border
  ) +
 # BS / BS nodes
  geom_point(
    data = bs_nodes,
    aes(
      x = x,
      y = y
    ),
    shape = 21,
    size = 6.0,
    stroke = 0.45,
    fill = color_bs,
    color = color_border
  ) +
 # RS / RS labels
  geom_text(
    data = rs_nodes,
    aes(
      x = x - 0.13,
      y = y,
      label = Label
    ),
    hjust = 1,
    size = 3.5,
    family = "serif",
    color = color_border
  ) +
 # ,
 # Taxon labels in italics
  geom_text(
    data = taxa_nodes,
    aes(
      x = x,
      y = y - 0.22,
      label = Label
    ),
    hjust = 0.5,
    vjust = 1,
    size = 3.25,
    family = "serif",
    fontface = "italic",
    color = color_border
  ) +
 # BS / BS labels
  geom_text(
    data = bs_nodes,
    aes(
      x = x + 0.13,
      y = y,
      label = Label
    ),
    hjust = 0,
    size = 3.5,
    family = "serif",
    color = color_border
  ) +
 # / Column titles
  annotate(
    geom = "text",
    x = 1,
    y = 5.48,
    label = "Rhizosphere soil (RS)",
    size = 4.0,
    family = "serif",
    fontface = "bold",
    color = color_border
  ) +
  annotate(
    geom = "text",
    x = 2.5,
    y = 5.48,
    label = "Yield-associated taxa",
    size = 4.0,
    family = "serif",
    fontface = "bold",
    color = color_border
  ) +
  annotate(
    geom = "text",
    x = 4,
    y = 5.48,
    label = "Bulk soil (BS)",
    size = 4.0,
    family = "serif",
    fontface = "bold",
    color = color_border
  ) +
 # / Ordering note
  annotate(
    geom = "text",
    x = 1,
    y = 0.67,
    label = "Higher contribution  →  lower contribution",
    size = 2.9,
    family = "serif",
    fontface = "italic",
    color = "#555555"
  ) +
  annotate(
    geom = "text",
    x = 4,
    y = 0.67,
    label = "Higher contribution  →  lower contribution",
    size = 2.9,
    family = "serif",
    fontface = "italic",
    color = "#555555"
  ) +
  coord_cartesian(
    xlim = c(
      0.25,
      4.75
    ),
    ylim = c(
      0.25,
      6.00
    ),
    clip = "off"
  ) +
  labs(
    caption = paste0(
      "Only global BH-adjusted FDR < ",
      fdr_threshold,
      " associations are displayed."
    )
  ) +
  theme_void(
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
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8.5
    ),
    plot.caption = element_text(
      hjust = 0.5,
      size = 8.5,
      color = "#444444",
      margin = margin(
        t = 8
      )
    ),
    plot.margin = margin(
      t = 10,
      r = 18,
      b = 10,
      l = 18
    )
  )
# ============================================================
# 16. PDF,SVGPNG
# Export PDF, SVG and PNG
# ============================================================
pdf_file <- file.path(
  output_dir,
  "RS_four_taxa_BS_network.pdf"
)
svg_file <- file.path(
  output_dir,
  "RS_four_taxa_BS_network.svg"
)
png_file <- file.path(
  output_dir,
  "RS_four_taxa_BS_network.png"
)
ggsave(
  filename = pdf_file,
  plot = network_plot,
  width = 11.2,
  height = 6.6,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "transparent"
)
grDevices::svg(
  filename = svg_file,
  width = 11.2,
  height = 6.6,
  bg = "transparent",
  onefile = TRUE
)
print(
  network_plot
)
grDevices::dev.off()
ggsave(
  filename = png_file,
  plot = network_plot,
  width = 11.2,
  height = 6.6,
  units = "in",
  dpi = 600,
  bg = "white"
)
# ============================================================
# 17.
# Export node table
# ============================================================
write_csv(
  nodes,
  file.path(
    output_dir,
    "Network_nodes.csv"
  )
)
# ============================================================
# 18.
# Completion messages
# ============================================================
message(
  "============================================================"
)
message(
  "."
)
message(
  "RS:",
  paste(
    rs_indicator_order,
    collapse = ";"
  )
)
message(
  "BS:",
  paste(
    bs_indicator_order,
    collapse = ";"
  )
)
message(
  "4:",
  paste(
    focal_taxa,
    collapse = ";"
  )
)
message(
  ":",
  nrow(edges)
)
message(
  "PDF:",
  normalizePath(pdf_file)
)
message(
  "SVG:",
  normalizePath(svg_file)
)
message(
  "PNG:",
  normalizePath(png_file)
)
message(
  ":",
  normalizePath(
    file.path(
      output_dir,
      "RS_BS_four_taxa_Spearman_all_results.csv"
    )
  )
)
message(
  "============================================================"
)
