# ======================================================================
# FIGURE S1
# Panel A: NMDS trajectories for BS, RS and RT compartments
# Panel B: ternary enrichment profiles for Control and Managed treatments
# This script is the author-supplied workflow corresponding to Figure S1.
# Legacy internal labels referring to Figure 2 were renamed to Figure S1
# for consistency with the final manuscript.
# ======================================================================
# Sanitized upload-ready script
# Figure S1B: ternary enrichment profiles
rm(list = ls())
suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(ggplot2)
  library(ggtern)
  library(grid)
})
# === 2. CSV Genus/Taxonomy ===
genus_file <- "genus_abundance.csv"
if (!file.exists(genus_file)) stop("Input file not found: ", genus_file)
data <- read.table(genus_file, header = TRUE, sep = ",", check.names = FALSE)
rownames(data) <- data$Genus
data <- data[, -1]
# === 3. ===
check_cols <- function(data, prefix) {
  pat <- function(suf) paste0("^", prefix, "[-_]*", suf)
  bs_cols <- grep(pat("BS"), colnames(data), value = TRUE)
  rs_cols <- grep(pat("RS"), colnames(data), value = TRUE)
  rt_cols <- grep(pat("RT"), colnames(data), value = TRUE)
  cat("====", prefix, "====\n",
      "BS cols (", length(bs_cols), "): ", paste(bs_cols, collapse=", "), "\n",
      "RS cols (", length(rs_cols), "): ", paste(rs_cols, collapse=", "), "\n",
      "RT cols (", length(rt_cols), "): ", paste(rt_cols, collapse=", "), "\n", sep = "")
  invisible(list(BS=bs_cols, RS=rs_cols, RT=rt_cols))
}
# === 4. BS/RS/RT ===
classify_enrichment <- function(df, fold = 1.2, min_prop = 0.4, delta = 0.05) {
  vals <- as.matrix(df[, c("BS","RS","RT")])
  labs <- c("BS","RS","RT")
  cls <- apply(vals, 1, function(x){
    x[is.na(x)] <- 0
    if (sum(x) <= 0) return("None")
    ord <- order(x, decreasing = TRUE)
    max_ix <- ord[1]; snd_ix <- ord[2]
    cond_fold  <- x[max_ix] >= x[snd_ix] * fold
    cond_prop  <- x[max_ix] >= min_prop
    cond_delta <- (x[max_ix] - x[snd_ix]) >= delta
    if (cond_fold && cond_prop && cond_delta) paste0(labs[max_ix], "_enriched") else "None"
  })
  df$Enrichment <- factor(cls, levels = c("BS_enriched","RS_enriched","RT_enriched","None"))
  df
}
# === 5. -> -> ===
extract_avg <- function(data, prefix, fold = 1.2, min_prop = 0.4, delta = 0.05) {
  pick <- function(pat, label) {
    idx <- grep(pat, colnames(data))
    if (length(idx) == 0) {
      stop(sprintf("%s : pattern = '%s'\n '%sBS'/'%sRS'/'%sRT' '-'  '_' ",
                   label, pat, prefix, prefix, prefix))
    }
    list(vals = rowMeans(data[, idx, drop = FALSE], na.rm = TRUE),
         cols = colnames(data)[idx], n = length(idx))
  }
  BS <- pick(paste0("^", prefix, "[-_]*BS"), paste0(prefix, "-BS"))
  RS <- pick(paste0("^", prefix, "[-_]*RS"), paste0(prefix, "-RS"))
  RT <- pick(paste0("^", prefix, "[-_]*RT"), paste0(prefix, "-RT"))
  df <- data.frame(BS_raw = BS$vals, RS_raw = RS$vals, RT_raw = RT$vals,
                   Abundance = rowMeans(cbind(BS$vals, RS$vals, RT$vals), na.rm = TRUE))
  rownames(df) <- rownames(data)
  prop <- t(apply(df[, c("BS_raw","RS_raw","RT_raw")], 1, function(x){
    x[is.na(x)] <- 0; s <- sum(x); if (s == 0) c(0,0,0) else x/s
  }))
  colnames(prop) <- c("BS","RS","RT")
  df <- cbind(df, prop)
  df <- tibble::rownames_to_column(df, var = "Taxonomy")
  df$RowID <- seq_len(nrow(df))
  df$Label <- df$Taxonomy
  attr(df, "match_info") <- list(
    BS_cols = BS$cols, RS_cols = RS$cols, RT_cols = RT$cols,
    n_BS = BS$n, n_RS = RS$n, n_RT = RT$n
  )
  classify_enrichment(df, fold = fold, min_prop = min_prop, delta = delta)
}
# === 6. near / outside ===
draw_enrichment_plot <- function(
    df, group_label,
    top_n = 10,
    mode = c("outside","near"),
 # outside
    pad = 0.08, gap = 0.045,              #
    min_gap = 0.030,                      #
    lanes = 1, lane_pad = 0.045,          #
    side_cap = Inf,                       # Inf
    elbow_len = 0.14, elbow_t = 0.65,     #
 # near
    label_offset = 0.10, jitter_sd = 0.008,
    min_dist = 0.030, step_inc = 0.010, max_iter = 80,
    tick_len = 0.25, show_tick = TRUE,
    line_color = "#E67E22",               # +
    text_size = 3.0
){
  mode <- match.arg(mode)
 # BS=RS=RT==
  enrich_colors <- c(
    "BS_enriched"="#7A0019",
    "RS_enriched"="#E67E22",
    "RT_enriched"="#8DAA1A",
    "None"       ="grey80"
  )
 # top_n
  pool <- df[df$Enrichment != "None",
             c("RowID","Label","Taxonomy","BS","RS","RT","Abundance","Enrichment")]
  pool <- pool[order(pool$Abundance, decreasing = TRUE), ]
  ann  <- if (nrow(pool)) pool[seq_len(min(top_n, nrow(pool))), ] else pool
  message(sprintf(
    "[%s]  = %d %d %s",
    group_label, nrow(pool), nrow(ann),
    if (nrow(ann)) paste(ann$Taxonomy, collapse = ", ") else ""
  ))
 # ===== OUTSIDE + =====
  if (mode == "outside" && nrow(ann) > 0) {
 # >= min_gap
    pack_1d <- function(x, min_gap, lo = 1e-6, hi = 1-1e-6){
      if(length(x) <= 1) return(pmin(pmax(x, lo), hi))
      x <- sort(x); x[1] <- max(x[1], lo)
      for(i in 2:length(x)) x[i] <- max(x[i], x[i-1] + min_gap)
      over <- x[length(x)] - hi
      if (over > 0) x <- x - over
      pmin(pmax(x, lo), hi)
    }
    place_side <- function(dat, side){
      if(nrow(dat)==0) return(dat[0,])
      if(side=="BS")      ord <- order(dat$RS)
      else if(side=="RS") ord <- order(dat$RT)
      else                ord <- order(dat$BS)
      dat <- dat[ord, , drop=FALSE]
      if (is.finite(side_cap) && nrow(dat) > side_cap)
        dat <- dat[seq_len(side_cap), , drop = FALSE]
      dat$lane <- (seq_len(nrow(dat)) - 1) %% lanes
      base <- 1 - 1e-6
      if(side=="BS"){
        dat$BS_lab <- 1 + pad + dat$lane * lane_pad
        raw <- seq(0, 1, length.out = nrow(dat)+2); raw <- raw[-c(1,length(raw))]
        raw <- pack_1d(raw, min_gap)
        dat$RS_lab <- pmin(base, pmax(1e-6, raw))
        dat$RT_lab <- 1 - dat$RS_lab
      } else if(side=="RS"){
        dat$RS_lab <- 1 + pad + dat$lane * lane_pad
        raw <- seq(0, 1, length.out = nrow(dat)+2); raw <- raw[-c(1,length(raw))]
        raw <- pack_1d(raw, min_gap)
        dat$RT_lab <- pmin(base, pmax(1e-6, raw))
        dat$BS_lab <- 1 - dat$RT_lab
      } else { # RT
        dat$RT_lab <- 1 + pad + dat$lane * lane_pad
        raw <- seq(0, 1, length.out = nrow(dat)+2); raw <- raw[-c(1,length(raw))]
        raw <- pack_1d(raw, min_gap)
        dat$BS_lab <- pmin(base, pmax(1e-6, raw))
        dat$RS_lab <- 1 - dat$BS_lab
      }
      dat
    }
    ann_bs <- place_side(subset(ann, Enrichment=="BS_enriched"), "BS")
    ann_rs <- place_side(subset(ann, Enrichment=="RS_enriched"), "RS")
    ann_rt <- place_side(subset(ann, Enrichment=="RT_enriched"), "RT")
    ann_out <- rbind(ann_bs, ann_rs, ann_rt)
    ann_out$BS_mid <- NA; ann_out$RS_mid <- NA; ann_out$RT_mid <- NA
    for (i in seq_len(nrow(ann_out))) {
      side <- if (ann_out$Enrichment[i]=="BS_enriched") "BS"
      else if (ann_out$Enrichment[i]=="RS_enriched") "RS" else "RT"
      if (side == "BS") {
        bs1 <- ann_out$BS[i] + elbow_len
        rs1 <- ann_out$RS[i] - elbow_len/2
        rt1 <- ann_out$RT[i] - elbow_len/2
      } else if (side == "RS") {
        bs1 <- ann_out$BS[i] - elbow_len/2
        rs1 <- ann_out$RS[i] + elbow_len
        rt1 <- ann_out$RT[i] - elbow_len/2
      } else { # RT
        bs1 <- ann_out$BS[i] - elbow_len/2
        rs1 <- ann_out$RS[i] - elbow_len/2
        rt1 <- ann_out$RT[i] + elbow_len
      }
      ann_out$BS_mid[i] <- bs1 + (ann_out$BS_lab[i] - bs1) * elbow_t
      ann_out$RS_mid[i] <- rs1 + (ann_out$RS_lab[i] - rs1) * elbow_t
      ann_out$RT_mid[i] <- rt1 + (ann_out$RT_lab[i] - rt1) * elbow_t
    }
  }
 # ===== NEAR =====
  if (mode == "near" && nrow(ann) > 0) {
    Cb <- 1/3; Cr <- 1/3; Ct <- 1/3
    ann$vBS <- ann$BS - Cb; ann$vRS <- ann$RS - Cr; ann$vRT <- ann$RT - Ct
    ann$vn  <- sqrt(ann$vBS^2 + ann$vRS^2 + ann$vRT^2) + 1e-9
    ann$uBS <- ann$vBS/ann$vn; ann$uRS <- ann$vRS/ann$vn; ann$uRT <- ann$vRT/ann$vn
    ann$BS_lab <- ann$BS + label_offset*ann$uBS + rnorm(nrow(ann),0,jitter_sd)
    ann$RS_lab <- ann$RS + label_offset*ann$uRS + rnorm(nrow(ann),0,jitter_sd)
    ann$RT_lab <- ann$RT + label_offset*ann$uRT + rnorm(nrow(ann),0,jitter_sd)
    clamp01 <- function(x) pmin(pmax(x, 1e-5), 1-1e-5)
    placed <- list()
    for (i in seq_len(nrow(ann))) {
      bs <- ann$BS_lab[i]; rs <- ann$RS_lab[i]; rt <- ann$RT_lab[i]; it <- 0
      repeat {
        too_close <- FALSE
        if (length(placed)) {
          d <- sapply(placed, function(p) sqrt((bs-p[1])^2 + (rs-p[2])^2 + (rt-p[3])^2))
          if (any(d < min_dist)) too_close <- TRUE
        }
        if (!too_close || it >= max_iter) break
        bs <- clamp01(bs + step_inc*ann$uBS[i])
        rs <- clamp01(rs + step_inc*ann$uRS[i])
        rt <- clamp01(1 - bs - rs)
        it <- it + 1
      }
      ann$BS_lab[i] <- bs; ann$RS_lab[i] <- rs; ann$RT_lab[i] <- rt
      placed[[length(placed)+1]] <- c(bs, rs, rt)
    }
    if (show_tick) {
      ann$BS_tick <- ann$BS + tick_len*(ann$BS_lab - ann$BS)
      ann$RS_tick <- ann$RS + tick_len*(ann$RS_lab - ann$RS)
      ann$RT_tick <- ann$RT + tick_len*(ann$RT_lab - ann$RT)
    }
  }
  p <- ggtern(data = df, aes(x = RS, y = BS, z = RT)) +
    geom_mask() +
    geom_point(aes(fill = Enrichment, size = Abundance),
               shape = 21, color = "black", stroke = 0.25, alpha = 0.85) +
    scale_fill_manual(values = enrich_colors, drop = FALSE,
                      name = "Enrichment",
                      labels = c("BS_enriched"="Enriched in BS",
                                 "RS_enriched"="Enriched in RS",
                                 "RT_enriched"="Enriched in RT",
                                 "None"="Unenriched")) +
    scale_size(range = c(1.8, 7.2), name = "Abundance (%)") +
    labs(T = "BS", L = "RS", R = "RT") +
    theme_bw(base_size = 12) + theme_showarrows() +
    theme(
      legend.position  = "right",
      panel.grid.major = element_line(size = 0.25, linetype = "dashed", color = "grey75"),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(size = 12, face = "bold"),
      plot.margin      = margin(6, 28, 6, 12)
    )
 # ===== & =====
  if (nrow(ann) > 0) {
    if (mode == "outside") {
      p <- p +
        geom_point(
          data = ann,
          aes(x = RS, y = BS, z = RT, size = Abundance),
          inherit.aes = FALSE, shape = 21, fill = NA,
          color = line_color, stroke = 0.9
        ) +
        geom_segment(
          data = ann_out,
          aes(x = RS, y = BS, z = RT,
              xend = RS_mid, yend = BS_mid, zend = RT_mid),
          inherit.aes = FALSE, color = line_color,
          linewidth = 0.5, lineend = "round"
        ) +
        geom_segment(
          data = ann_out,
          aes(x = RS_mid, y = BS_mid, z = RT_mid,
              xend = RS_lab, yend = BS_lab, zend = RT_lab),
          inherit.aes = FALSE, color = line_color,
          linewidth = 0.5, lineend = "round"
        ) +
        geom_text(
          data = ann_out,
          aes(x = RS_lab, y = BS_lab, z = RT_lab, label = Label),
          inherit.aes = FALSE, color = line_color,
          size = text_size, fontface = "bold",
          position = ggplot2::position_identity()  #   ggtern
        )
    } else { # near
      p <- p +
        geom_point(
          data = ann,
          aes(x = RS, y = BS, z = RT, size = Abundance),
          inherit.aes = FALSE, shape = 21, fill = NA,
          color = line_color, stroke = 0.9
        ) +
        { if (show_tick) list(
          geom_segment(
            data = ann,
            aes(x = RS, y = BS, z = RT,
                xend = RS_tick, yend = BS_tick, zend = RT_tick),
            inherit.aes = FALSE, color = line_color,
            linewidth = 0.45, lineend = "round"
          )
        ) } +
        geom_text(
          data = ann,
          aes(x = RS_lab, y = BS_lab, z = RT_lab, label = Label),
          inherit.aes = FALSE, color = line_color,
          size = text_size, fontface = "bold",
          position = ggplot2::position_identity()  #
        )
    }
  }
  p
}
# === 7. panel ===
to_noclip <- function(p) {
  gt <- ggplotGrob(p)
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  gt
}
# === 8. ===
export_enriched_detail <- function(df, group_label, outfile) {
  mi <- attr(df, "match_info"); stopifnot(!is.null(mi))
  sub <- subset(df, Enrichment != "None",
                select = c("Taxonomy","BS_raw","RS_raw","RT_raw","BS","RS","RT","Abundance","Enrichment"))
  sub$n_BS <- mi$n_BS; sub$n_RS <- mi$n_RS; sub$n_RT <- mi$n_RT
  sub <- sub[order(sub$Enrichment, -sub$Abundance), ]
  write.csv(sub, outfile, row.names = FALSE)
  message(sprintf("[%s]  : %s (=%d)", group_label, outfile, nrow(sub)))
  sub
}
# === 9. ===
check_cols(data, "G"); check_cols(data, "H")
data_g <- extract_avg(data, "G", fold = 1.2, min_prop = 0.4, delta = 0.05)
data_h <- extract_avg(data, "H", fold = 1.2, min_prop = 0.4, delta = 0.05)
# 10
plot_g <- draw_enrichment_plot(
  data_g, "G_group",
  top_n = 10, mode = "outside",
  pad = 0.08, gap = 0.045, min_gap = 0.030,
  lanes = 2, lane_pad = 0.045, side_cap = Inf,
  elbow_len = 0.14, elbow_t = 0.65,
  text_size = 3.0
)
plot_h <- draw_enrichment_plot(
  data_h, "H_group",
  top_n = 10, mode = "outside",
  pad = 0.08, gap = 0.045, min_gap = 0.030,
  lanes = 2, lane_pad = 0.045, side_cap = Inf,
  elbow_len = 0.14, elbow_t = 0.65,
  text_size = 3.0
)
grid::grid.newpage(); grid::grid.draw(to_noclip(plot_g))
grid::grid.newpage(); grid::grid.draw(to_noclip(plot_h))
ggsave("FigS1B_Control_ternary_enrichment.pdf", to_noclip(plot_g), width = 10, height = 9, dpi = 600)
ggsave("FigS1B_Managed_ternary_enrichment.pdf", to_noclip(plot_h), width = 10, height = 9, dpi = 600)
export_enriched_detail(data_g, "G_group", "G_enriched_detail_check.csv")
export_enriched_detail(data_h, "H_group", "H_enriched_detail_check.csv")
write.csv(data_g, "data_g_full.csv", row.names = FALSE)
write.csv(data_h, "data_h_full.csv", row.names = FALSE)
# Figure S1A: NMDS trajectories
# NMDS X/Y +
suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ggforce)
  library(ggnewscale)
  library(magrittr)
})
nmds_file <- "nmds_coordinates.csv"
if (!file.exists(nmds_file)) stop("Input file not found: ", nmds_file)
df0 <- readr::read_csv(nmds_file, show_col_types = FALSE)
stopifnot(all(c("sample","NMDS1","NMDS2","Group") %in% names(df0)))
# ====== Group ======
parsed <- str_match(df0$Group, "^([GH])(BS|RS|RT)[-_]?(\\d+|CK)$")
df <- df0 %>%
  mutate(
    Treatment = parsed[,2],
    MainGroup = parsed[,3],
    Time_raw  = parsed[,4]
  )
bad <- which(is.na(df$Treatment) | is.na(df$MainGroup) | is.na(df$Time_raw))
if (length(bad) > 0) {
  message("  Group ",
          paste(unique(df$Group[bad]), collapse = ", "))
  df <- df[-bad, , drop = FALSE]
}
df <- df %>%
  mutate(
    Time = case_when(
      is.na(Time_raw) ~ NA_real_,
      toupper(Time_raw) == "CK" ~ 0,
      grepl("^\\d+$", Time_raw) ~ as.numeric(Time_raw),
      TRUE ~ NA_real_
    ),
    Series = paste0(Treatment, MainGroup)
  )
# ====== / ======
time_for_scale <- ifelse(is.na(df$Time), 0, df$Time)
df <- df %>%
  mutate(
    size       = scales::rescale(time_for_scale, to = c(3, 9),
                                 from = range(time_for_scale, na.rm = TRUE)),
    alpha_fill = scales::rescale(time_for_scale, to = c(0.35, 0.75),
                                 from = range(time_for_scale, na.rm = TRUE))
  )
col_treat <- c(G = "#6984af", H = "#cf6248")  #  vs
# ====== / ======
compute_curve_df <- function(d) {
  dff <- d %>%
    group_by(Series, Time) %>%
    summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop") %>%
    arrange(Series, Time) %>%
    group_by(Series) %>%
    slice(c(1, n())) %>%
    mutate(pos = c("start","end")) %>%
    ungroup()
  start_df <- dff %>% filter(pos == "start") %>%
    select(Series, NMDS1_start = NMDS1, NMDS2_start = NMDS2)
  end_df <- dff %>% filter(pos == "end") %>%
    select(Series, NMDS1_end = NMDS1, NMDS2_end = NMDS2)
  curve_df <- left_join(start_df, end_df, by = "Series") %>%
    filter(!(NMDS1_start == NMDS1_end & NMDS2_start == NMDS2_end))
  series_levels <- unique(curve_df$Series)
  base_curves <- c(0.3, -0.3, 0.5, -0.5, 0.7, -0.7, 0.4, -0.4, 0.6, -0.6)
  curve_map <- tibble(
    Series = series_levels,
    curvature = rep(base_curves, length.out = length(series_levels))
  )
  curve_df %>%
    left_join(curve_map, by = "Series") %>%
    mutate(curvature = ifelse(is.na(curvature), 0.4, curvature))
}
compute_equal_limits <- function(d, pad = 0.10) {
  rx <- range(d$NMDS1, na.rm = TRUE)
  ry <- range(d$NMDS2, na.rm = TRUE)
  cx <- mean(rx); cy <- mean(ry)
  sx <- diff(rx); sy <- diff(ry)
  if (!is.finite(sx) || sx == 0) sx <- 1e-6
  if (!is.finite(sy) || sy == 0) sy <- 1e-6
  span <- max(sx, sy) * (1 + pad*2)
  list(
    xlim = c(cx - span/2, cx + span/2),
    ylim = c(cy - span/2, cy + span/2)
  )
}
make_plot_one <- function(main_group, data_all, outfile = NULL,
                          xlim = NULL, ylim = NULL, ratio = 1, auto_pad = 0.10) {
  d <- data_all %>% filter(MainGroup == main_group)
  if (nrow(d) == 0) return(invisible(NULL))
  curve_df2 <- compute_curve_df(d)
  lims <- compute_equal_limits(d, pad = auto_pad)
  if (is.null(xlim)) xlim <- lims$xlim
  if (is.null(ylim)) ylim <- lims$ylim
  curve_layers <- purrr::pmap(
    curve_df2,
    function(Series, NMDS1_start, NMDS2_start, NMDS1_end, NMDS2_end, curvature) {
      Treat <- substr(Series, 1, 1)
      ggplot2::geom_curve(
        aes(x = NMDS1_start, y = NMDS2_start,
            xend = NMDS1_end, yend = NMDS2_end),
        data = tibble(Series, NMDS1_start, NMDS2_start, NMDS1_end, NMDS2_end),
        color = col_treat[[Treat]],
        curvature = curvature,
        arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
        linewidth = 0.9,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    }
  )
  p <- ggplot(d, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(fill = Treatment, size = size, alpha = alpha_fill),
               shape = 21, color = "grey25", stroke = 0.3) +
    scale_size_identity() +
    scale_alpha_identity() +
    scale_fill_manual(values = col_treat) +
    guides(size = "none", alpha = "none",
           fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    theme_test() +
    theme(
      legend.title = element_blank(),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = paste0("NMDS Trajectory  ", main_group)) +
    coord_fixed(ratio = ratio, xlim = xlim, ylim = ylim) +
    curve_layers
  if (!is.null(outfile)) {
    ggsave(outfile, p, width = 8, height = 6, dpi = 300)
    message(" ", outfile)
  }
  return(p)
}
# Output files are written to the working directory.
p_bs <- make_plot_one("BS", df, "FigS1A_NMDS_Trajectory_BS.pdf",
                      xlim = c(-0.015, -0.005), ylim = c(-0.005, 0.005), ratio = 1)
p_rs <- make_plot_one("RS", df, "FigS1A_NMDS_Trajectory_RS.pdf",
                      xlim = c(-0.015, -0.005), ylim = c(-0.003, 0.007), ratio = 1)
p_rt <- make_plot_one("RT", df, "FigS1A_NMDS_Trajectory_RT.pdf",
                      xlim = c(-0.1, 0.1), ylim = c(-0.12, 0.09), ratio = 1)
# ----------------------------------------------------------------------
# Final reproducibility record
# ----------------------------------------------------------------------
capture.output(
  sessionInfo(),
  file = "sessionInfo_FigS1.txt"
)
cat("\nFigure S1 analysis completed.\n")
