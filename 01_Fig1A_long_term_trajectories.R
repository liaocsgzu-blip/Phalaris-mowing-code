# ============================================================
# Figure 1A: long-term trajectories
# Upload-ready version with portable relative paths.
# Usage
# Rscript 01_Fig1A_long_term_trajectories.R [input.xlsx] [output_dir]
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else "main_plant_soil_data.xlsx"
output_dir <- if (length(args) >= 2) args[[2]] else "Fig1A_output"
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (dir.exists("work/Rlib")) {
  local_lib <- normalizePath("work/Rlib", winslash = "/", mustWork = TRUE)
  .libPaths(c(local_lib, .libPaths()))
}
suppressWarnings(suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
}))
palette <- c("Unmown" = "#55B6A8", "Mown" = "#E47B6B")
shape_values <- c("Unmown" = 16, "Mown" = 17)
line_values <- c("Unmown" = "22", "Mown" = "solid")
years <- c(0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16)
theme_nature <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      axis.title = element_text(size = 6.5, colour = "black"),
      axis.text = element_text(size = 5.8, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.key.width = unit(8, "mm"),
      legend.spacing.x = unit(1.2, "mm"),
      legend.position = "top",
      legend.box.margin = margin(0, 0, 0, 0),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.3, face = "bold"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 5.8, colour = "#444444"),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank(),
      panel.spacing = unit(2.5, "mm"),
      plot.margin = margin(2.5, 3.5, 2.5, 3.5)
    )
}
theme_set(theme_nature())
soil <- read_excel(input_file, sheet = "Soil")
plant <- read_excel(input_file, sheet = "Plant")
carbon <- read_excel(input_file, sheet = "CarbonUse")
names(soil) <- c(
  "soil_type", "treatment_code", "year", "total_n", "organic_matter",
  "total_p", "total_k", "alkali_n", "available_k", "available_p", "pH"
)
names(plant) <- c(
  "sample_id", "treatment_year", "year", "yield_t_ha", "root_mass_index",
  "root_activity", "soil_compaction"
)
names(carbon) <- c(
  "soil_type", "treatment_code", "year", "mbc", "mbn", "mbc_mbn",
  "growth_total", "growth_bacteria", "growth_fungi",
  "resp_total", "resp_bacteria", "resp_fungi",
  "cue", "bacterial_cue", "fungal_cue", "turnover_rate"
)
decode_treatment <- function(x) {
  factor(ifelse(x == "H", "Mown", "Unmown"), levels = c("Unmown", "Mown"))
}
soil <- soil |>
  mutate(
    treatment = decode_treatment(treatment_code),
    year = as.numeric(year)
  )
plant <- plant |>
  mutate(
    treatment_code = sub("[0-9].*$", "", treatment_year),
    treatment = decode_treatment(treatment_code),
    year = as.numeric(year)
  )
carbon <- carbon |>
  mutate(
    treatment = decode_treatment(treatment_code),
    year = as.numeric(year)
  )
summary_ci <- function(data, id_cols, variables) {
  data |>
    pivot_longer(all_of(variables), names_to = "variable", values_to = "value") |>
    group_by(across(all_of(c(id_cols, "treatment"))), variable) |>
    summarise(
      n = sum(!is.na(value)),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      se = sd / sqrt(n),
      ci95 = qt(0.975, df = n - 1) * se,
      .groups = "drop"
    )
}
yearwise_tests <- function(data, strata, variables) {
  data |>
    pivot_longer(all_of(variables), names_to = "variable", values_to = "value") |>
    group_by(across(all_of(strata)), variable) |>
    group_modify(~ {
      d <- .x |> filter(!is.na(value))
      means <- tapply(d$value, d$treatment, mean)
      tt <- t.test(value ~ treatment, data = d)
      tibble(
        estimate_mown_minus_unmown = unname(means[["Mown"]] - means[["Unmown"]]),
        p_value = tt$p.value
      )
    }) |>
    ungroup() |>
    group_by(variable) |>
    mutate(p_adj_bh = p.adjust(p_value, method = "BH")) |>
    ungroup()
}
plant_vars <- c("yield_t_ha", "root_mass_index", "root_activity", "soil_compaction")
soil_vars <- c("total_n", "organic_matter", "total_p", "total_k", "alkali_n", "available_k", "available_p", "pH")
carbon_vars <- c("cue", "bacterial_cue", "fungal_cue", "turnover_rate")
plant_summary <- summary_ci(plant, c("year"), plant_vars)
soil_summary <- summary_ci(soil, c("soil_type", "year"), soil_vars)
carbon_summary <- summary_ci(carbon, c("soil_type", "year"), carbon_vars)
plant_tests <- yearwise_tests(plant, c("year"), plant_vars)
soil_tests <- yearwise_tests(soil, c("soil_type", "year"), soil_vars)
carbon_tests <- yearwise_tests(carbon, c("soil_type", "year"), carbon_vars)
trajectory_plot <- function(
    raw_data, summary_data, test_data, variable_name, y_label, title,
    soil_compartment = NULL, facet_soil = FALSE, x_breaks = years,
    show_legend = TRUE, significance = TRUE
) {
  raw_long <- raw_data |>
    pivot_longer(all_of(variable_name), names_to = "variable", values_to = "value") |>
    filter(variable == variable_name)
  sum_df <- summary_data |> filter(variable == variable_name)
  test_df <- test_data |> filter(variable == variable_name)
  if (!is.null(soil_compartment)) {
    raw_long <- raw_long |> filter(soil_type == soil_compartment)
    sum_df <- sum_df |> filter(soil_type == soil_compartment)
    test_df <- test_df |> filter(soil_type == soil_compartment)
  }
  sig_df <- tibble()
  if (significance) {
    sig_df <- test_df |>
      filter(p_adj_bh < 0.05) |>
      left_join(
        sum_df |>
          group_by(across(any_of(c("soil_type", "year")))) |>
          summarise(y_base = max(mean + ci95, na.rm = TRUE), .groups = "drop"),
        by = intersect(c("soil_type", "year"), names(test_df))
      ) |>
      group_by(across(any_of("soil_type"))) |>
      mutate(
        span = max(sum_df$mean + sum_df$ci95, na.rm = TRUE) -
          min(sum_df$mean - sum_df$ci95, na.rm = TRUE),
        y_sig = y_base + 0.075 * ifelse(span > 0, span, max(abs(sum_df$mean), na.rm = TRUE))
      ) |>
      ungroup()
  }
  p <- ggplot() +
    geom_ribbon(
      data = sum_df,
      aes(x = year, ymin = mean - ci95, ymax = mean + ci95, fill = treatment, group = treatment),
      alpha = 0.11, colour = NA
    ) +
    geom_point(
      data = raw_long,
      aes(x = year, y = value, colour = treatment, shape = treatment),
      alpha = 0.28, size = 0.8,
      position = position_jitter(width = 0.10, height = 0)
    ) +
    geom_line(
      data = sum_df,
      aes(x = year, y = mean, colour = treatment, linetype = treatment, group = treatment),
      linewidth = 0.72
    ) +
    geom_point(
      data = sum_df,
      aes(x = year, y = mean, colour = treatment, shape = treatment),
      size = 1.65, stroke = 0.3
    ) +
    scale_colour_manual(values = palette) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = shape_values) +
    scale_linetype_manual(values = line_values) +
    scale_x_continuous(breaks = x_breaks, minor_breaks = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
    labs(x = "Years since introduction", y = y_label, title = title) +
    theme_nature() +
    guides(fill = "none")
  if (nrow(sig_df) > 0) {
    p <- p + geom_text(
      data = sig_df,
      aes(x = year, y = y_sig, label = "*"),
      inherit.aes = FALSE, size = 2.1, family = "Arial", colour = "black", vjust = 0
    )
  }
  if (facet_soil) {
    p <- p + facet_wrap(~soil_type, nrow = 1, labeller = as_labeller(c(BS = "Bulk soil", RS = "Rhizosphere")))
  }
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}
save_pub_r <- function(plot, stem, width_mm, height_mm, dpi_tiff = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(
    file.path(output_dir, paste0(stem, ".svg")),
    width = width_in, height = height_in, bg = "white", system_fonts = list(sans = "Arial")
  )
  print(plot)
  dev.off()
  grDevices::cairo_pdf(
    file.path(output_dir, paste0(stem, ".pdf")),
    width = width_in, height = height_in, family = "Arial", bg = "white"
  )
  print(plot)
  dev.off()
  ragg::agg_tiff(
    file.path(output_dir, paste0(stem, ".tiff")),
    width = width_in, height = height_in, units = "in",
    res = dpi_tiff, background = "white", compression = "lzw"
  )
  print(plot)
  dev.off()
  ragg::agg_png(
    file.path(output_dir, paste0(stem, ".png")),
    width = width_in, height = height_in, units = "in",
    res = 300, background = "white"
  )
  print(plot)
  dev.off()
}
# Figure 1: plant performance.
p1a <- trajectory_plot(
  plant, plant_summary, plant_tests, "yield_t_ha",
  "Aboveground yield (t ha\u207b\u00b9)", "Yield", show_legend = TRUE
)
p1b <- trajectory_plot(
  plant, plant_summary, plant_tests, "root_mass_index",
  "Root-mass index", "Root mass", show_legend = FALSE
)
p1c <- trajectory_plot(
  plant, plant_summary, plant_tests, "root_activity",
  "Root activity (\u03bcg TTC g\u207b\u00b9 FW h\u207b\u00b9)", "Root activity", show_legend = FALSE
)
p1d <- trajectory_plot(
  plant, plant_summary, plant_tests, "soil_compaction",
  "Soil compaction (MPa)", "Soil compaction", show_legend = FALSE
)
fig1 <- p1a / (p1b | p1c | p1d) +
  plot_layout(heights = c(1.15, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top", plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(fig1, "Fig1_Plant_performance", 183, 125)
# Figure 2: soil organic matter and nitrogen.
p2a <- trajectory_plot(
  soil, soil_summary, soil_tests, "organic_matter",
  "Soil organic matter (g kg\u207b\u00b9)", "Rhizosphere", "RS", show_legend = TRUE
)
p2b <- trajectory_plot(
  soil, soil_summary, soil_tests, "organic_matter",
  "Soil organic matter (g kg\u207b\u00b9)", "Bulk soil", "BS", show_legend = FALSE
)
p2c <- trajectory_plot(
  soil, soil_summary, soil_tests, "total_n",
  "Total nitrogen (g kg\u207b\u00b9)", "Rhizosphere", "RS", show_legend = FALSE
)
p2d <- trajectory_plot(
  soil, soil_summary, soil_tests, "total_n",
  "Total nitrogen (g kg\u207b\u00b9)", "Bulk soil", "BS", show_legend = FALSE
)
fig2 <- (p2a | p2b) / (p2c | p2d) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top", plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(fig2, "Fig2_Soil_organic_matter_and_N", 183, 125)
# Figure 3: microbial carbon use.
short_years <- c(0, 4, 8, 12, 16)
p3a <- trajectory_plot(
  carbon, carbon_summary, carbon_tests, "cue",
  "Community CUE", "Community CUE", facet_soil = TRUE, x_breaks = short_years, show_legend = TRUE
)
p3b <- trajectory_plot(
  carbon, carbon_summary, carbon_tests, "bacterial_cue",
  "Bacterial CUE", "Bacterial CUE", facet_soil = TRUE, x_breaks = short_years, show_legend = FALSE
)
p3c <- trajectory_plot(
  carbon, carbon_summary, carbon_tests, "fungal_cue",
  "Fungal CUE", "Fungal CUE", facet_soil = TRUE, x_breaks = short_years, show_legend = FALSE
)
p3d <- trajectory_plot(
  carbon, carbon_summary, carbon_tests, "turnover_rate",
  "Turnover rate", "Microbial turnover", facet_soil = TRUE, x_breaks = short_years, show_legend = FALSE
)
fig3 <- (p3a | p3b) / (p3c | p3d) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top", plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(fig3, "Fig3_Microbial_CUE_and_turnover", 183, 135)
# Supplementary heatmap: log2 response ratio for the mown-unmown contrast.
soil_long <- soil |>
  pivot_longer(all_of(soil_vars), names_to = "variable", values_to = "value")
effect_data <- soil_long |>
  group_by(soil_type, year, variable) |>
  group_modify(~ {
    d <- .x
    x1 <- d$value[d$treatment == "Mown"]
    x0 <- d$value[d$treatment == "Unmown"]
    tt <- t.test(x1, x0)
    tibble(log2_response_ratio = log2(mean(x1) / mean(x0)), p_value = tt$p.value)
  }) |>
  ungroup() |>
  group_by(variable) |>
  mutate(p_adj_bh = p.adjust(p_value, method = "BH")) |>
  ungroup()
variable_labels <- c(
  organic_matter = "Organic matter",
  total_n = "Total N",
  total_p = "Total P",
  total_k = "Total K",
  alkali_n = "Alkali-hydrolysable N",
  available_k = "Available K",
  available_p = "Available P",
  pH = "pH"
)
effect_data <- effect_data |>
  mutate(
    variable_label = factor(
      variable_labels[variable],
      levels = rev(unname(variable_labels))
    ),
    sig = ifelse(p_adj_bh < 0.05, "*", "")
  )
effect_limit <- max(abs(effect_data$log2_response_ratio), na.rm = TRUE)
effect_data <- effect_data |>
  mutate(sig_colour = ifelse(abs(log2_response_ratio) > 0.55 * effect_limit, "white", "black"))
heatmap_panel <- function(compartment, title) {
  ggplot(effect_data |> filter(soil_type == compartment), aes(x = factor(year), y = variable_label, fill = log2_response_ratio)) +
    geom_tile(colour = "white", linewidth = 0.28) +
    geom_text(aes(label = sig, colour = sig_colour), family = "Arial", size = 2.1) +
    scale_fill_gradient2(
      low = "#3B6FB6", mid = "white", high = "#D55E00",
      midpoint = 0, limits = c(-effect_limit, effect_limit),
      name = "log\u2082 response ratio\nMown / unmown"
    ) +
    scale_colour_identity() +
    labs(x = "Years since introduction", y = NULL, title = title) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 0),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 6, colour = "black")
    )
}
figs1 <- heatmap_panel("RS", "Rhizosphere") | heatmap_panel("BS", "Bulk soil")
figs1 <- figs1 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "right", plot.tag = element_text(size = 8, face = "bold"))
# Additional diagnostic output; not part of the final Supplementary Figure S1.
save_pub_r(figs1, "Diagnostic_soil_nutrient_effects", 183, 115)
# Traceable source data.
write.csv(
  plant_summary |> left_join(plant_tests, by = c("year", "variable")),
  file.path(output_dir, "Fig1_source_data.csv"), row.names = FALSE
)
write.csv(
  soil_summary |>
    filter(variable %in% c("organic_matter", "total_n")) |>
    left_join(soil_tests, by = c("soil_type", "year", "variable")),
  file.path(output_dir, "Fig2_source_data.csv"), row.names = FALSE
)
write.csv(
  carbon_summary |> left_join(carbon_tests, by = c("soil_type", "year", "variable")),
  file.path(output_dir, "Fig3_source_data.csv"), row.names = FALSE
)
write.csv(effect_data, file.path(output_dir, "FigS1_source_data.csv"), row.names = FALSE)
notes <- c(
  "Treatment mapping used for all figures: H = Mown; G = Unmown.",
  "Biological replicates: n = 5 per treatment × year; soil and carbon panels also stratify by soil compartment.",
  "Points: raw observations. Lines/large symbols: group means. Shaded bands: 95% t confidence intervals.",
  "Asterisks: yearwise Welch t-test, Benjamini-Hochberg adjusted P < 0.05 within each variable.",
  "No observations were removed, imputed, or smoothed.",
  "The source workbook does not define the H/G mapping; confirm this before manuscript submission.",
  "The source root-mass header includes '100*'; the physical unit must be confirmed before final labelling.",
  "Soil organic matter is plotted as supplied and is not relabelled as soil organic carbon."
)
writeLines(notes, file.path(output_dir, "README_figure_notes.txt"), useBytes = TRUE)
if (!is.null(warnings())) {
  capture.output(warnings(), file = file.path(output_dir, "R_warnings.txt"))
} else {
  writeLines("No plotting warnings.", file.path(output_dir, "R_warnings.txt"))
}
cat("Generated figures in:", normalizePath(output_dir, winslash = "/"), "\n")
