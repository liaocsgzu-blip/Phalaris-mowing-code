args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  input <- args[[1]]
  output_dir <- args[[2]]
} else {
  input <- "main_plant_soil_data.xlsx"
  output_dir <- "Fig1BC_output"
  message("Using default relative paths.")
}
if (!file.exists(input)) {
  stop("Input file not found: ", input)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (dir.exists("work/Rlib")) {
  lib <- normalizePath("work/Rlib", winslash = "/", mustWork = TRUE)
  .libPaths(c(lib, .libPaths()))
}
required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "patchwork",
  "MuMIn", "relaimpo", "svglite", "ragg"
)
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them first, then rerun the script."
  )
}
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(MuMIn)
  library(relaimpo)
  library(svglite)
  library(ragg)
})
suppressPackageStartupMessages(library(dplyr))
soil <- read_excel(input, sheet = "Soil")
plant <- read_excel(input, sheet = "Plant")
carbon <- read_excel(input, sheet = "CarbonUse")
names(soil) <- c(
  "soil_type", "treatment", "year", "total_N", "organic_matter",
  "total_P", "total_K", "alkaline_N", "available_K", "available_P",
  "pH", "SOC"
)
names(plant) <- c(
  "sample_id", "treatment_year", "year", "yield", "root_mass",
  "root_activity", "soil_compaction"
)
names(carbon) <- c(
  "soil_type", "treatment", "year", "MBC", "MBN", "MBC_MBN",
  "growth_total", "growth_bacteria", "growth_fungi",
  "resp_total", "resp_bacteria", "resp_fungi",
  "CUE", "CUE_bacteria", "CUE_fungi", "turnover"
)
soil <- soil |>
  group_by(soil_type, treatment, year) |>
  mutate(rep = row_number()) |>
  ungroup()
carbon <- carbon |>
  group_by(soil_type, treatment, year) |>
  mutate(rep = row_number()) |>
  ungroup()
plant <- plant |>
  mutate(
    treatment = substr(treatment_year, 1, 1),
    rep = as.integer(sub(".*-", "", sample_id))
  )
dat <- soil |>
  inner_join(carbon, by = c("soil_type", "treatment", "year", "rep")) |>
  inner_join(
    plant |> dplyr::select(treatment, year, rep, yield, root_mass, root_activity, soil_compaction),
    by = c("treatment", "year", "rep")
  ) |>
  arrange(soil_type, treatment, year, rep)
if (nrow(dat) != 220L || anyNA(dat)) {
  stop("Unexpected merged data structure: expected 220 complete rows.")
}
if (any(duplicated(dat[c("soil_type", "treatment", "year", "rep")]))) {
  stop("Merged data contain duplicate observation keys.")
}
predictors_all <- c(
  "total_N", "organic_matter", "total_P", "total_K", "alkaline_N",
  "available_K", "available_P", "pH", "SOC",
  "MBC", "MBN", "MBC_MBN", "growth_total", "growth_bacteria",
  "resp_total", "resp_bacteria",
  "CUE", "CUE_bacteria", "turnover"
)
variable_labels <- c(
  total_N = "TN", organic_matter = "OM", total_P = "TP", total_K = "TK",
  alkaline_N = "AN", available_K = "AK", available_P = "AP", pH = "pH",
  SOC = "SOC", MBC = "MBC", MBN = "MBN", MBC_MBN = "MBC:MBN",
  growth_total = "Growth", growth_bacteria = "Bact. growth",
  growth_fungi = "Fungal growth", resp_total = "Respiration",
  resp_bacteria = "Bact. resp.", resp_fungi = "Fungal resp.",
  CUE = "CUE", CUE_bacteria = "Bact. CUE", CUE_fungi = "Fungal CUE",
  turnover = "Turnover", root_mass = "Root mass",
  root_activity = "Root activity", soil_compaction = "Compaction"
)
variable_dictionary <- data.frame(
  variable = names(variable_labels),
  label = unname(variable_labels),
  description = c(
    "Total nitrogen", "Organic matter", "Total phosphorus", "Total potassium",
    "Alkaline-hydrolysable nitrogen", "Available potassium", "Available phosphorus",
    "Soil pH", "Soil organic carbon", "Microbial biomass carbon",
    "Microbial biomass nitrogen", "Microbial biomass C:N ratio",
    "Total microbial growth rate", "Bacterial growth rate", "Fungal growth rate",
    "Total microbial respiration rate", "Bacterial respiration rate",
    "Fungal respiration rate", "Microbial carbon-use efficiency",
    "Bacterial carbon-use efficiency", "Fungal carbon-use efficiency",
    "Microbial turnover rate", "Root mass", "Root activity", "Soil compaction"
  ),
  stringsAsFactors = FALSE
) |>
  filter(variable %in% predictors_all)
calc_vif <- function(df, vars) {
  setNames(vapply(vars, function(v) {
    others <- setdiff(vars, v)
    fit <- lm(reformulate(others, response = v), data = df)
    r2 <- summary(fit)$r.squared
    if (!is.finite(r2) || r2 >= 1 - 1e-12) Inf else 1 / (1 - r2)
  }, numeric(1)), vars)
}
# Use one common predictor set so RS and BS importance values remain comparable.
current <- predictors_all
vif_audit <- list()
removed <- data.frame(
  variable = character(), removal_step = integer(), max_vif = numeric(),
  rs_vif = numeric(), bs_vif = numeric()
)
step <- 0L
repeat {
  step <- step + 1L
  z_rs <- as.data.frame(scale(dat |> filter(soil_type == "RS") |> dplyr::select(all_of(current))))
  z_bs <- as.data.frame(scale(dat |> filter(soil_type == "BS") |> dplyr::select(all_of(current))))
  v_rs <- calc_vif(z_rs, current)
  v_bs <- calc_vif(z_bs, current)
  vmax <- pmax(v_rs, v_bs)
  vif_audit[[step]] <- bind_rows(
    data.frame(step = step, soil_type = "RS", variable = names(v_rs), vif = as.numeric(v_rs)),
    data.frame(step = step, soil_type = "BS", variable = names(v_bs), vif = as.numeric(v_bs))
  )
  if (max(vmax) <= 10 || length(current) <= 2) break
  drop_var <- names(which.max(vmax))
  removed <- bind_rows(
    removed,
    data.frame(
      variable = drop_var, removal_step = step, max_vif = vmax[[drop_var]],
      rs_vif = v_rs[[drop_var]], bs_vif = v_bs[[drop_var]]
    )
  )
  current <- setdiff(current, drop_var)
}
predictors <- current
vif_final <- bind_rows(
  data.frame(soil_type = "RS", variable = names(v_rs), vif = as.numeric(v_rs)),
  data.frame(soil_type = "BS", variable = names(v_bs), vif = as.numeric(v_bs))
)
screening <- variable_dictionary |>
  left_join(removed, by = "variable") |>
  mutate(
    status = if_else(variable %in% predictors, "Retained (final VIF <= 10)", "Excluded during iterative VIF screening"),
    final_max_vif = pmax(
      unname(v_rs[match(variable, names(v_rs))]),
      unname(v_bs[match(variable, names(v_bs))])
    )
  )
options(na.action = "na.fail")
analyse_stratum <- function(st) {
  raw <- dat |> filter(soil_type == st) |> dplyr::select(yield, all_of(predictors))
  z <- as.data.frame(scale(raw))
  global <- lm(reformulate(predictors, response = "yield"), data = z)
  dd <- dredge(global, rank = "BIC", trace = FALSE)
  top_index <- which(dd$delta < 4)
  top_models <- get.models(dd, subset = delta < 4)
  top_table <- as.data.frame(dd[top_index, , drop = FALSE])
  top_table$model_rank <- seq_len(nrow(top_table))
  top_table$model_id <- paste0(st, "_M", top_table$model_rank)
  top_table$weight_top <- top_table$weight / sum(top_table$weight)
  top_table$soil_type <- st
  inclusion <- bind_rows(lapply(seq_along(top_models), function(i) {
    terms_i <- attr(terms(top_models[[i]]), "term.labels")
    data.frame(
      soil_type = st,
      model_rank = i,
      model_id = paste0(st, "_M", i),
      variable = predictors,
      included = predictors %in% terms_i
    )
  }))
  model_metrics <- top_table |>
    transmute(
      soil_type, model_rank, model_id, df, logLik, BIC, delta,
      weight_global = weight, weight_top
    )
  lmg <- calc.relimp(global, type = "lmg", rela = FALSE)$lmg
  coefs <- coef(global)[predictors]
  top_presence <- inclusion |>
    left_join(model_metrics |> dplyr::select(model_id, weight_top), by = "model_id") |>
    group_by(variable) |>
    summarise(top_model_weight = sum(weight_top[included]), .groups = "drop")
  importance <- data.frame(
    soil_type = st,
    variable = predictors,
    lmg_percent = 100 * as.numeric(lmg[predictors]),
    standardized_beta = as.numeric(coefs[predictors])
  ) |>
    left_join(top_presence, by = "variable")
  fit_stats <- data.frame(
    soil_type = st,
    n = nrow(z),
    predictors_retained = length(predictors),
    global_R2 = summary(global)$r.squared,
    global_adjusted_R2 = summary(global)$adj.r.squared,
    global_BIC = BIC(global),
    models_evaluated = nrow(dd),
    supported_models_delta_BIC_lt_4 = nrow(top_table),
    best_model_BIC = min(dd$BIC)
  )
  list(
    standardized_data = z,
    dredge = dd,
    top_table = top_table,
    inclusion = inclusion,
    model_metrics = model_metrics,
    importance = importance,
    fit_stats = fit_stats
  )
}
results <- lapply(c("RS", "BS"), analyse_stratum)
names(results) <- c("RS", "BS")
importance <- bind_rows(lapply(results, `[[`, "importance")) |>
  left_join(variable_dictionary, by = "variable")
inclusion <- bind_rows(lapply(results, `[[`, "inclusion")) |>
  left_join(variable_dictionary |> dplyr::select(variable, label), by = "variable")
model_metrics <- bind_rows(lapply(results, `[[`, "model_metrics"))
fit_stats <- bind_rows(lapply(results, `[[`, "fit_stats"))
analysis_data <- dat |>
  dplyr::select(soil_type, treatment, year, rep, yield, all_of(predictors_all))
write.csv(analysis_data, file.path(output_dir, "matched_analysis_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bind_rows(vif_audit), file.path(output_dir, "vif_audit_all_steps.csv"), row.names = FALSE)
write.csv(vif_final, file.path(output_dir, "vif_final.csv"), row.names = FALSE)
write.csv(screening, file.path(output_dir, "variable_screening.csv"), row.names = FALSE)
write.csv(variable_dictionary, file.path(output_dir, "variable_dictionary.csv"), row.names = FALSE)
write.csv(importance, file.path(output_dir, "relative_importance.csv"), row.names = FALSE)
write.csv(fit_stats, file.path(output_dir, "model_fit_summary.csv"), row.names = FALSE)
for (st in c("RS", "BS")) {
  write.csv(
    results[[st]]$model_metrics,
    file.path(output_dir, paste0("supported_models_", st, ".csv")),
    row.names = FALSE
  )
  write.csv(
    as.data.frame(results[[st]]$dredge),
    file.path(output_dir, paste0("all_models_", st, ".csv")),
    row.names = FALSE
  )
}
plot_labels <- variable_labels[predictors]
palette <- c(
  "#4E79A7", "#59A14F", "#F28E2B", "#E15759", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6", "#6F8FC9",
  "#A0CBE8", "#8CD17D", "#D4A6C8", "#C49A6C", "#7BA7A5", "#E3B46F",
  "#8B8FB9", "#A6A6A6"
)[seq_along(predictors)]
names(palette) <- predictors
metrics_for_labels <- model_metrics |>
  mutate(soil_type = factor(soil_type, levels = c("RS", "BS"))) |>
  arrange(soil_type, model_rank) |>
  mutate(
    model_label = sprintf(
      "M%d   BIC %.1f   Δ %.2f   w %.2f",
      model_rank, BIC, delta, weight_top
    )
  )
display_n <- 12L
supported_counts <- setNames(fit_stats$supported_models_delta_BIC_lt_4, fit_stats$soil_type)
shown_counts <- setNames(
  pmin(display_n, as.numeric(supported_counts)),
  names(supported_counts)
)
strip_labels <- c(
  RS = sprintf("RS | all %d supported models shown", shown_counts[["RS"]]),
  BS = sprintf("BS | all %d supported models shown", shown_counts[["BS"]])
)
inclusion_plot <- inclusion |>
  filter(model_rank <= display_n) |>
  left_join(metrics_for_labels |> dplyr::select(soil_type, model_rank, model_label),
            by = c("soil_type", "model_rank")) |>
  mutate(
    soil_type = factor(soil_type, levels = c("RS", "BS")),
    label = factor(label, levels = unname(plot_labels)),
    model_label = factor(model_label, levels = rev(unique(model_label))),
    fill_key = if_else(included, variable, "Absent")
  )
fill_values <- c(Absent = "#F2F2F2", palette)
theme_pub <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 7.2, face = "bold"),
    plot.tag = element_text(size = 9, face = "bold"),
    panel.grid = element_blank()
  )
p_a <- ggplot(inclusion_plot, aes(x = label, y = model_label, fill = fill_key)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_manual(values = fill_values, guide = "none") +
  facet_grid(soil_type ~ ., scales = "free_y", space = "free_y",
             labeller = as_labeller(strip_labels)) +
  scale_x_discrete(position = "top", drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 48, hjust = 0, vjust = 0, size = 6.2),
    axis.text.y = element_text(size = 5.8),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.spacing.y = unit(4, "pt"),
    plot.margin = margin(3, 4, 3, 3)
  )
importance_order <- importance |>
  group_by(variable, label) |>
  summarise(mean_importance = mean(lmg_percent), .groups = "drop") |>
  arrange(mean_importance) |>
  pull(label)
importance_plot <- importance |>
  mutate(
    soil_type = factor(soil_type, levels = c("RS", "BS")),
    label = factor(label, levels = importance_order)
  )
p_b <- ggplot(importance_plot, aes(x = lmg_percent, y = label, fill = variable)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.1f", lmg_percent)),
    hjust = -0.12, size = 2.05, family = "Arial"
  ) +
  scale_fill_manual(values = palette, guide = "none") +
  facet_wrap(
    ~soil_type, nrow = 1,
    labeller = as_labeller(c(RS = "Rhizosphere soil (RS)", BS = "Bulk soil (BS)"))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15)),
    breaks = scales::pretty_breaks(5)
  ) +
  labs(x = "Variance explained by each predictor (%)", y = NULL) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 6.1),
    axis.title.x = element_text(size = 7),
    panel.spacing.x = unit(12, "pt"),
    plot.margin = margin(3, 4, 3, 3)
  )
fig <- p_a / p_b +
  plot_layout(heights = c(1.0, 1.45)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "Arial"))
base <- file.path(output_dir, "yield_predictor_importance")
width_mm <- 183
height_mm <- 185
w <- width_mm / 25.4
h <- height_mm / 25.4
svglite::svglite(paste0(base, ".svg"), width = w, height = h, system_fonts = list(Arial = "Arial"))
print(fig)
dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = "Arial")
print(fig)
dev.off()
ragg::agg_tiff(
  paste0(base, ".tiff"), width = w, height = h, units = "in",
  res = 600, background = "white"
)
print(fig)
dev.off()
ragg::agg_png(
  paste0(base, ".png"), width = w, height = h, units = "in",
  res = 300, background = "white"
)
print(fig)
dev.off()
make_split_figure <- function(st) {
  st_title <- if (st == "RS") "Rhizosphere soil (RS)" else "Bulk soil (BS)"
  ref_label_map <- c(
    total_N = "TN", organic_matter = "OM", total_P = "TP",
    total_K = "TK", alkaline_N = "AN", available_K = "AK",
    available_P = "AP", pH = "pH", SOC = "SOC",
    MBC = "MBC", MBN = "MBN", MBC_MBN = "MBC:MBN",
    growth_total = "Growth", growth_bacteria = "BG",
    resp_total = "Resp.", resp_bacteria = "BR",
    CUE = "CUE", CUE_bacteria = "BCUE", turnover = "Turnover"
  )
  inc_st <- inclusion_plot |>
    filter(soil_type == st) |>
    mutate(
      label_ref = factor(
        unname(ref_label_map[variable]),
        levels = unname(ref_label_map[predictors])
      ),
      model_rank_factor = factor(
        model_rank,
        levels = rev(sort(unique(model_rank)))
      )
    ) |>
    droplevels()
  imp_st <- importance_plot |>
    filter(soil_type == st) |>
    droplevels()
  p_matrix <- ggplot(inc_st, aes(x = label_ref, y = model_rank_factor, fill = fill_key)) +
    geom_tile(colour = "white", linewidth = 0.28) +
    scale_fill_manual(values = fill_values, guide = "none") +
    scale_x_discrete(position = "top", drop = FALSE) +
    labs(tag = "A", x = NULL, y = st) +
    theme_pub +
    theme(
      plot.tag = element_text(size = 9, face = "bold"),
      plot.tag.position = c(0, 1.12),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 5.7),
      axis.text.y = element_blank(),
      axis.title.y = element_text(size = 6.2, angle = 0, vjust = 0.5),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.margin = margin(8, 1, 2, 5)
    )
  metrics_st <- metrics_for_labels |>
    filter(as.character(soil_type) == st) |>
    mutate(
      model_rank_factor = factor(
        model_rank,
        levels = rev(sort(unique(model_rank)))
      ),
      df_txt = sprintf("%d", df),
      logLik_txt = sprintf("%.2f", logLik),
      BIC_txt = sprintf("%.1f", BIC),
      delta_txt = sprintf("%.2f", delta),
      weight_txt = sprintf("%.2f", weight_top)
    ) |>
    dplyr::select(model_rank_factor, df_txt, logLik_txt, BIC_txt, delta_txt, weight_txt) |>
    pivot_longer(
      -model_rank_factor,
      names_to = "metric",
      values_to = "value"
    ) |>
    mutate(
      metric = factor(
        metric,
        levels = c("df_txt", "logLik_txt", "BIC_txt", "delta_txt", "weight_txt"),
        labels = c("df", "LogLik", "BIC", "Delta", "Weight")
      )
    )
  p_table <- ggplot(metrics_st, aes(x = metric, y = model_rank_factor, label = value)) +
    geom_text(size = 2.05, family = "Arial") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_void(base_family = "Arial") +
    theme(
      axis.text.x = element_text(size = 5.5, colour = "black"),
      plot.margin = margin(8, 4, 2, 1)
    )
  p_top <- (p_matrix | p_table) +
    plot_layout(widths = c(4.4, 1.9))
  imp_stack <- imp_st |>
    mutate(order_id = match(variable, predictors)) |>
    arrange(order_id) |>
    mutate(
      xmin = lag(cumsum(lmg_percent), default = 0),
      xmax = cumsum(lmg_percent),
      xmid = (xmin + xmax) / 2,
      label_short = unname(ref_label_map[variable]),
      narrow = lmg_percent < 3.0,
      narrow_index = cumsum(narrow),
      label_y = if_else(
        narrow,
        if_else(narrow_index %% 2 == 0, 0.79, 0.21),
        0.50
      )
    )
  p_importance <- ggplot(imp_stack) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = 0.34, ymax = 0.66, fill = variable),
      colour = "white", linewidth = 0.28
    ) +
    scale_fill_manual(values = palette, guide = "none") +
    geom_segment(
      data = imp_stack |> filter(narrow),
      aes(
        x = xmid, xend = xmid,
        y = if_else(label_y > 0.5, 0.66, 0.34),
        yend = if_else(label_y > 0.5, 0.73, 0.27)
      ),
      linewidth = 0.25, colour = "#555555"
    ) +
    geom_text(
      aes(x = xmid, y = label_y, label = label_short),
      size = 2.0, family = "Arial"
    ) +
    scale_x_continuous(
      limits = c(0, 80),
      expand = c(0, 0),
      breaks = seq(0, 80, 20)
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(tag = "B", x = "Variance explained in the model (%)", y = NULL) +
    theme_pub +
    theme(
      plot.tag = element_text(size = 9, face = "bold"),
      plot.tag.position = c(0, 1.02),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_line(linewidth = 0.35),
      axis.title.x = element_text(size = 7),
      plot.margin = margin(5, 10, 4, 5)
    )
  split_fig <- (p_top / p_importance) +
    plot_layout(heights = c(1.15, 1.0)) +
    plot_annotation(
      title = st_title,
      theme = theme(
        plot.title = element_text(
          family = "Arial", size = 8, face = "bold",
          hjust = 0, margin = margin(b = 2)
        )
      )
    )
  split_base <- file.path(output_dir, paste0("yield_reference_style_importance_", st))
  split_width_mm <- 183
  split_height_mm <- 125
  split_w <- split_width_mm / 25.4
  split_h <- split_height_mm / 25.4
  svglite::svglite(paste0(split_base, ".svg"), width = split_w, height = split_h,
                   system_fonts = list(Arial = "Arial"))
  print(split_fig)
  dev.off()
  grDevices::cairo_pdf(paste0(split_base, ".pdf"), width = split_w, height = split_h,
                       family = "Arial")
  print(split_fig)
  dev.off()
  ragg::agg_tiff(paste0(split_base, ".tiff"), width = split_w, height = split_h,
                 units = "in", res = 600, background = "white")
  print(split_fig)
  dev.off()
  ragg::agg_png(paste0(split_base, ".png"), width = split_w, height = split_h,
                units = "in", res = 300, background = "white")
  print(split_fig)
  dev.off()
}
invisible(lapply(c("RS", "BS"), make_split_figure))
report <- c(
  "Yield predictor importance analysis",
  "===================================",
  "",
  "Core conclusion: Predictor importance is evaluated separately for rhizosphere (RS) and bulk soil (BS) to avoid duplicating each plant-yield observation.",
  "Response: yield (t/ha).",
  "Candidate predictors: measured soil and microbial carbon-use indicators after excluding all fungi-specific variables; plant traits, soil compaction, treatment and year are not ranked as predictors.",
  "Standardization: response and retained predictors z-standardized within each soil stratum.",
  "Collinearity: a common predictor set was obtained by iterative VIF screening; at each step the variable with the largest VIF across RS and BS was removed until both strata had VIF <= 10.",
  "Model selection: all subsets of retained predictors; BIC ranking; supported set defined as delta BIC < 4.",
  "Relative importance: relaimpo::calc.relimp(type = 'lmg'); values are percentage points of total response variance and sum to global R-squared x 100.",
  "",
  paste0("Retained predictors (n = ", length(predictors), "): ", paste(predictors, collapse = ", ")),
  "",
  capture.output(print(fit_stats, row.names = FALSE)),
  "",
  "Reviewer-risk notes:",
  "- Results are associational, not causal.",
  "- RS and BS are analysed separately because they share the same plant-yield observations.",
  "- VIF screening removes redundant indicators; excluded variables remain listed in variable_screening.csv.",
  "- Relative importance can be sensitive to the candidate predictor set and strong ecological gradients."
)
writeLines(report, file.path(output_dir, "analysis_report.txt"), useBytes = TRUE)
message("Analysis complete: ", normalizePath(output_dir, winslash = "/"))
