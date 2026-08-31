# ============================================================
# RT + PLS:R
# RT metatranscriptomic functional co-expression and cross-year PLS
# 1. CSV39RT
# 2. RT(SEM)
# 3. /CLR
# 4. ,()
# 5. Functional co-expression network construction
# 6. ,PLS
# 7. kMEVIP,Nature.
# - "",
# - WGCNA
# - NA0.NA, na_as_zero <- FALSE.
# ============================================================
options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260727)
# ---------------------------
# 0.
# ---------------------------
input_dir  <- "."                  # CSV
output_dir <- file.path(input_dir, "RT_R_resultsMR")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
single_panel_dir <- file.path(output_dir, "Single_panels")
dir.create(single_panel_dir, showWarnings = FALSE, recursive = TRUE)
supplementary_panel_dir <- file.path(output_dir, "Supplementary_single_panels")
dir.create(supplementary_panel_dir, showWarnings = FALSE, recursive = TRUE)
na_as_zero <- TRUE
prevalence_cutoff <- 0.30
top_n_network_features <- 250
soft_power_override <- 10          # ;NULL
candidate_powers <- 1:12
min_mean_connectivity <- 5
tree_cut_height <- 0.905
min_module_size <- 12
n_preservation_permutations <- 300
max_pls_components <- 4
# Input-file discovery
find_input_file <- function(pattern, label) {
  x <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(x) == 0) stop("Missing input file: ", label)
  sort(x)[1]
}
file_meta    <- find_input_file("^sample_mapping.*\\.csv$", "sample_mapping*.csv")
file_species <- find_input_file("^rt_species_abundance.*\\.csv$", "rt_species_abundance*.csv")
file_plant   <- find_input_file("^plant_traits_yield.*\\.csv$", "plant_traits_yield*.csv")
file_metab   <- find_input_file("^root_metabolome_summary.*\\.csv$", "root_metabolome_summary*.csv")
file_mg      <- find_input_file("^metagenome_selected_functions.*\\.csv$", "metagenome_selected_functions*.csv")
file_mt      <- find_input_file("^metatranscriptome_selected_functions.*\\.csv$", "metatranscriptome_selected_functions*.csv")
# ---------------------------
# 1. R
# ---------------------------
cran_packages <- c(
  "data.table", "stringr", "ggplot2", "patchwork", "ggrepel",
  "pls", "sandwich", "lmtest", "scales", "RColorBrewer",
  "ggdendro", "svglite"
)
installed <- rownames(installed.packages())
missing_packages <- setdiff(cran_packages, installed)
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org", dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(pls)
  library(sandwich)
  library(lmtest)
  library(scales)
  library(RColorBrewer)
  library(ggdendro)
  library(svglite)
})
theme_set(
  theme_classic(base_size = 9) +
    theme(
      text = element_text(family = "sans", colour = "black"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      plot.subtitle = element_text(size = 8),
      plot.tag = element_text(face = "bold", size = 11),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.45, "cm")
    )
)
# ---------------------------
# 2.
# ---------------------------
z_score <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
# Population-SD standardization used for the SEM-defined composite states.
# SEM:.
z_pop <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  s <- sqrt(mean((x - mu)^2, na.rm = TRUE))
  if (!is.finite(s) || s < 1e-10) return(rep(0, length(x)))
  (x - mu) / s
}
rank_z <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  z_score(r)
}
safe_bh <- function(p) {
  out <- rep(NA_real_, length(p))
  keep <- is.finite(p)
  if (any(keep)) out[keep] <- p.adjust(p[keep], method = "BH")
  out
}
r2_score <- function(observed, predicted) {
  1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
}
safe_scale_train_test <- function(x_train, x_test) {
  mu <- colMeans(x_train, na.rm = TRUE)
  sigma <- apply(x_train, 2, sd, na.rm = TRUE)
  sigma[!is.finite(sigma) | sigma < 1e-10] <- 1
  train_scaled <- sweep(sweep(x_train, 2, mu, "-"), 2, sigma, "/")
  test_scaled  <- sweep(sweep(x_test,  2, mu, "-"), 2, sigma, "/")
  list(train = train_scaled, test = test_scaled, center = mu, scale = sigma)
}
robust_rank_association <- function(x, y, cell_factor) {
  dat <- data.frame(
    y = rank_z(y),
    x = rank_z(x),
    cell = factor(cell_factor)
  )
  fit <- lm(y ~ x + cell, data = dat)
  ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  beta <- unname(ct["x", "Estimate"])
  se <- unname(ct["x", "Std. Error"])
  p <- unname(ct["x", "Pr(>|t|)"])
  c(beta = beta, se = se, ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se, p = p)
}
robust_interaction_wald <- function(fit, pattern = ":") {
  coef_names <- names(coef(fit))
  idx <- grep(pattern, coef_names, fixed = TRUE)
  if (length(idx) == 0) return(c(statistic = NA, df = 0, p = NA))
  b <- coef(fit)[idx]
  V <- sandwich::vcovHC(fit, type = "HC3")[idx, idx, drop = FALSE]
  W <- as.numeric(t(b) %*% solve(V, b))
  c(statistic = W, df = length(idx), p = pchisq(W, df = length(idx), lower.tail = FALSE))
}
# ---------------------------
# 3.
# ---------------------------
meta <- fread(file_meta, encoding = "UTF-8", check.names = FALSE)
stopifnot("GeneSampleID" %in% names(meta))
sample_ids <- meta$GeneSampleID
species <- fread(file_species, encoding = "UTF-8", check.names = FALSE)
plant <- fread(file_plant, encoding = "UTF-8", check.names = FALSE)
metab <- fread(file_metab, encoding = "UTF-8", check.names = FALSE)
if (!all(sample_ids %in% names(species))) stop("16SRT.")
if (!all(sample_ids %in% names(metab))) stop("RT.")
# 3.1 RT:SEM
high_taxa <- c("Staphylococcus capitis", "Pseudomonas syncyanea")
low_taxa  <- c("Luteitalea pratensis", "Paraflavitalea soli")
if (!all(c(high_taxa, low_taxa) %in% species$Species)) {
  stop("16S4,.")
}
taxa_mat <- as.matrix(
  species[Species %in% c(high_taxa, low_taxa), ..sample_ids]
)
rownames(taxa_mat) <- species[Species %in% c(high_taxa, low_taxa), Species]
taxa_mat <- taxa_mat[c(high_taxa, low_taxa), , drop = FALSE]
storage.mode(taxa_mat) <- "double"
taxa_positive <- taxa_mat[is.finite(taxa_mat) & taxa_mat > 0]
taxa_pc <- min(taxa_positive, na.rm = TRUE) / 2
log_taxa <- log(taxa_mat + taxa_pc)
bacterial_state <- colMeans(log_taxa[high_taxa, , drop = FALSE]) -
  colMeans(log_taxa[low_taxa, , drop = FALSE])
# 3.2
a_plant <- match(sample_ids, plant$GeneSampleID)
if (anyNA(a_plant)) stop("RT.")
plant_sub <- plant[a_plant]
yield_value <- plant_sub[["Yield_t_ha"]]
# 3.3 :SEM
# 8 -> log1p -> SD -> -> PCA PC1
positive_metabolites <- c(
  "Com_1057_pos", # Sakuranetin
  "Com_659_pos",  # Hypoxanthine
  "Com_384_pos",  # Asp-Glu
  "Com_745_pos"   # L-Kynurenine
)
negative_metabolites <- c(
  "Com_701_neg",  # Uridine
  "Com_212_neg",  # Catalpol
  "Com_697_neg",  # Tyrosol
  "Com_146_neg"   # 8-O-Acetylharpagide
)
all_metabolites <- c(positive_metabolites, negative_metabolites)
if (!all(all_metabolites %in% metab$FeatureID)) {
  stop("8FeatureID.")
}
met_sub <- metab[FeatureID %in% all_metabolites]
met_mat <- as.matrix(met_sub[, ..sample_ids])
rownames(met_mat) <- met_sub$FeatureID
met_mat <- met_mat[all_metabolites, , drop = FALSE]
storage.mode(met_mat) <- "double"
# samples x metabolites
met_log_sample <- t(log1p(met_mat))
rownames(met_log_sample) <- sample_ids
colnames(met_log_sample) <- all_metabolites
met_z <- apply(met_log_sample, 2, z_pop)
met_z <- as.matrix(met_z)
rownames(met_z) <- sample_ids
colnames(met_z) <- all_metabolites
# Align the four negative metabolites so that higher values consistently
# represent the same metabolic direction used in the SEM.
met_z[, negative_metabolites] <- -met_z[, negative_metabolites, drop = FALSE]
met_pca <- prcomp(met_z, center = FALSE, scale. = FALSE)
metabolic_state <- as.numeric(met_pca$x[, 1])
# Orient PC1 toward a higher coordinated positive metabolic state.
if (cor(metabolic_state, rowMeans(met_z), use = "complete.obs") < 0) {
  metabolic_state <- -metabolic_state
  met_pca$rotation[, 1] <- -met_pca$rotation[, 1]
}
metabolic_state_variance <- met_pca$sdev[1]^2 / sum(met_pca$sdev^2)
metabolic_state_loadings <- data.frame(
  FeatureID = rownames(met_pca$rotation),
  PC1_loading = met_pca$rotation[, 1],
  stringsAsFactors = FALSE
)
traits <- data.frame(
  GeneSampleID = sample_ids,
  Yield = yield_value,
  Bacterial_state = bacterial_state,
  Metabolic_state = metabolic_state,
  TreatmentCode = meta$TreatmentCode,
  Year = meta$Year,
  Replicate = meta$Replicate,
  stringsAsFactors = FALSE
)
traits$Cell <- interaction(traits$TreatmentCode, traits$Year, drop = TRUE)
# ---------------------------
# 4. CLR
# ---------------------------
make_functional_unit <- function(database, feature_id, annotation) {
  out <- as.character(feature_id)
  idx <- database == "KEGG"
  if (any(idx)) {
    parts <- stringr::str_split_fixed(annotation[idx], " \\| ", 3)
    gene <- stringr::str_trim(stringr::str_extract(parts[, 2], "^[^;]+"))
    gene[is.na(gene)] <- ""
    out[idx] <- ifelse(gene == "", feature_id[idx], paste0(feature_id[idx], "|", gene))
  }
  idx <- database == "CAZy"
  if (any(idx)) out[idx] <- feature_id[idx]
  idx <- database %in% c("Pcyc", "Polyphenol")
  if (any(idx)) {
    parts <- stringr::str_split_fixed(annotation[idx], " \\| ", 6)
    candidate <- paste(parts[, 3], parts[, 4], parts[, 5], sep = "|")
    bad <- apply(parts[, 3:5, drop = FALSE], 1, function(z) any(is.na(z) | z == ""))
    candidate[bad] <- annotation[idx][bad]
    out[idx] <- candidate
  }
  idx <- database == "TCDB"
  if (any(idx)) {
    ann <- annotation[idx]
    gene <- stringr::str_match(ann, "GN=([A-Za-z0-9_-]+)")[, 2]
    tc <- stringr::str_match(ann, "gnl\\|TC-DB\\|[^|]+\\|([^|]+)")[, 2]
    candidate <- ifelse(
      !is.na(gene) & gene != "",
      paste0(gene, "|", ifelse(is.na(tc), feature_id[idx], tc)),
      ifelse(!is.na(tc) & tc != "", tc, feature_id[idx])
    )
    out[idx] <- candidate
  }
  out
}
aggregate_functional_units <- function(file_path, sample_ids, na_as_zero = TRUE) {
  message(":", basename(file_path))
  usecols <- c("Database", "Module", "FeatureID", "Annotation", sample_ids)
  dat <- fread(file_path, select = usecols, encoding = "UTF-8", check.names = FALSE)
  missing_sample_cols <- setdiff(sample_ids, names(dat))
  if (length(missing_sample_cols) > 0) {
    stop(":", paste(missing_sample_cols, collapse = ", "))
  }
  if (!na_as_zero && anyNA(dat[, ..sample_ids])) {
    stop("NA,na_as_zero=FALSE.NA.")
  }
  if (na_as_zero) {
    for (nm in sample_ids) set(dat, which(is.na(dat[[nm]])), nm, 0)
  }
  dat[, Unit := make_functional_unit(Database, FeatureID, Annotation)]
  agg <- dat[
    ,
    lapply(.SD, sum, na.rm = TRUE),
    by = .(Database, FunctionalModule = Module, Unit),
    .SDcols = sample_ids
  ]
  clr_blocks <- list()
  info_blocks <- list()
  for (db in unique(agg$Database)) {
    block <- agg[Database == db]
    mat <- as.matrix(block[, ..sample_ids])
    storage.mode(mat) <- "double"
    totals <- colSums(mat)
    totals[totals == 0] <- NA_real_
    relative <- sweep(mat, 2, totals, "/")
    relative[!is.finite(relative)] <- 0
    prevalence <- rowMeans(relative > 0)
    keep <- prevalence >= prevalence_cutoff
    block <- block[keep]
    relative <- relative[keep, , drop = FALSE]
    positive <- relative[relative > 0]
    pc <- if (length(positive) > 0) as.numeric(quantile(positive, 0.001, names = FALSE)) / 2 else 1e-10
    log_relative <- log(relative + pc)
    clr <- sweep(log_relative, 2, colMeans(log_relative), "-")
    key <- paste(block$Database, block$FunctionalModule, block$Unit, sep = "|||")
    rownames(clr) <- key
    clr_blocks[[db]] <- clr
    info_blocks[[db]] <- data.frame(
      FeatureKey = key,
      Database = block$Database,
      FunctionalModule = block$FunctionalModule,
      Unit = block$Unit,
      Prevalence = prevalence[keep],
      stringsAsFactors = FALSE
    )
  }
  clr_all <- do.call(rbind, clr_blocks)
  info_all <- do.call(rbind, info_blocks)
  variance <- apply(clr_all, 1, var)
  keep <- is.finite(variance) & variance > 1e-8
  list(
    matrix = clr_all[keep, , drop = FALSE],
    info = info_all[match(rownames(clr_all)[keep], info_all$FeatureKey), , drop = FALSE]
  )
}
mt_obj <- aggregate_functional_units(file_mt, sample_ids, na_as_zero)
mg_obj <- aggregate_functional_units(file_mg, sample_ids, na_as_zero)
mt_all <- mt_obj$matrix
mg_all <- mg_obj$matrix
feature_info <- mt_obj$info
rownames(feature_info) <- feature_info$FeatureKey
# ---------------------------
# 5.
# ---------------------------
mad_values <- apply(mt_all, 1, mad, constant = 1.4826)
network_keys <- names(sort(mad_values, decreasing = TRUE))[seq_len(min(top_n_network_features, length(mad_values)))]
X_network <- t(mt_all[network_keys, , drop = FALSE])
X_network <- apply(X_network, 2, z_score)
X_network <- as.matrix(X_network)
rownames(X_network) <- sample_ids
colnames(X_network) <- network_keys
cor_mat <- cor(X_network, use = "pairwise.complete.obs", method = "pearson")
diag(cor_mat) <- 0
scale_free_fit <- function(adjacency) {
  k <- rowSums(adjacency)
  h <- hist(k, breaks = max(10, floor(sqrt(length(k)))), plot = FALSE)
  keep <- h$counts > 0 & h$mids > 0
  if (sum(keep) < 5) return(c(R2 = NA, mean_k = mean(k)))
  fit <- lm(log10(h$counts[keep] / sum(h$counts)) ~ log10(h$mids[keep]))
  c(R2 = summary(fit)$r.squared, mean_k = mean(k))
}
power_table <- rbindlist(lapply(candidate_powers, function(power) {
  adj <- ((1 + cor_mat) / 2)^power
  diag(adj) <- 0
  fit <- scale_free_fit(adj)
  data.table(Power = power, ScaleFree_R2 = fit["R2"], MeanConnectivity = fit["mean_k"])
}))
if (is.null(soft_power_override)) {
  candidates <- power_table[MeanConnectivity >= min_mean_connectivity]
  soft_power <- candidates[order(-ScaleFree_R2, Power)][1, Power]
} else {
  soft_power <- soft_power_override
}
adjacency <- ((1 + cor_mat) / 2)^soft_power
adjacency[!is.finite(adjacency)] <- 0
adjacency <- (adjacency + t(adjacency)) / 2
diag(adjacency) <- 0
connectivity <- rowSums(adjacency)
shared_neighbors <- adjacency %*% adjacency
# TOMmin(k_i, k_j)
min_connectivity <- outer(
  connectivity,
  connectivity,
  FUN = function(x, y) pmin(x, y)
)
denominator <- min_connectivity + 1 - adjacency
denominator[!is.finite(denominator) | denominator <= 0] <- 1
TOM <- (shared_neighbors + adjacency) / denominator
TOM[!is.finite(TOM)] <- 0
TOM <- (TOM + t(TOM)) / 2
diag(TOM) <- 1
# pmax(0, pmin(1, 1 - TOM))
# pmin/pmax;dim,
# as.dist()/.
dissTOM <- 1 - TOM
dissTOM[!is.finite(dissTOM)] <- 1
dissTOM[dissTOM < 0] <- 0
dissTOM[dissTOM > 1] <- 1
dissTOM <- (dissTOM + t(dissTOM)) / 2
diag(dissTOM) <- 0
if (!is.matrix(dissTOM) || nrow(dissTOM) != ncol(dissTOM)) {
  stop(
    "dissTOM:",
    paste(dim(dissTOM), collapse = " × "),
    ";X_networkcor_mat."
  )
}
if (any(!is.finite(dissTOM))) {
  stop("dissTOMNA/NaN/Inf.")
}
message(
  ":X_network=", nrow(X_network), "×", ncol(X_network),
  ";cor_mat=", nrow(cor_mat), "×", ncol(cor_mat),
  ";dissTOM=", nrow(dissTOM), "×", ncol(dissTOM)
)
tree <- hclust(as.dist(dissTOM), method = "average")
raw_labels <- cutree(tree, h = tree_cut_height)
raw_sizes <- table(raw_labels)
clean_labels <- ifelse(raw_sizes[as.character(raw_labels)] >= min_module_size, raw_labels, 0)
valid_labels <- setdiff(unique(clean_labels), 0)
valid_labels <- valid_labels[order(sapply(valid_labels, function(z) sum(clean_labels == z)), decreasing = TRUE)]
label_map <- setNames(paste0("M", seq_along(valid_labels)), valid_labels)
module_assignment <- ifelse(clean_labels == 0, "Grey", label_map[as.character(clean_labels)])
names(module_assignment) <- colnames(X_network)
module_names <- paste0("M", seq_along(valid_labels))
base_module_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442")
module_colors <- setNames(base_module_colors[seq_along(module_names)], module_names)
module_colors <- c(module_colors, Grey = "#BDBDBD")
# ---------------------------
# 6. ()kME
# ---------------------------
eigengenes <- matrix(NA_real_, nrow = length(sample_ids), ncol = length(module_names),
                     dimnames = list(sample_ids, module_names))
loading_list <- list()
module_summary <- list()
member_results <- list()
for (module in module_names) {
  members <- names(module_assignment)[module_assignment == module]
  mat <- X_network[, members, drop = FALSE]
  pca <- prcomp(mat, center = FALSE, scale. = FALSE)
  me <- pca$x[, 1]
  loadings <- pca$rotation[, 1]
  member_cor <- apply(mat, 2, function(v) cor(v, me))
  if (mean(member_cor, na.rm = TRUE) < 0) {
    me <- -me
    loadings <- -loadings
    member_cor <- -member_cor
  }
  eigengenes[, module] <- me
  loading_list[[module]] <- loadings
  module_summary[[module]] <- data.frame(
    Module = module,
    N_features = length(members),
    PC1_variance_explained = pca$sdev[1]^2 / sum(pca$sdev^2),
    Mean_member_correlation = mean(member_cor, na.rm = TRUE)
  )
  member_p <- sapply(member_cor, function(r) {
    if (!is.finite(r) || abs(r) >= 1) return(0)
    tval <- r * sqrt((nrow(mat) - 2) / (1 - r^2))
    2 * pt(abs(tval), df = nrow(mat) - 2, lower.tail = FALSE)
  })
  member_q <- p.adjust(member_p, method = "BH")
  info <- feature_info[members, c("FeatureKey", "Database", "FunctionalModule", "Unit", "Prevalence"), drop = FALSE]
  info$Module <- module
  info$kME <- member_cor[members]
  info$kME_p <- member_p[members]
  info$kME_q <- member_q[members]
  member_results[[module]] <- info
}
module_summary <- rbindlist(module_summary)
module_members <- rbindlist(member_results)
grey_keys <- names(module_assignment)[module_assignment == "Grey"]
if (length(grey_keys) > 0) {
  grey <- feature_info[grey_keys, c("FeatureKey", "Database", "FunctionalModule", "Unit", "Prevalence"), drop = FALSE]
  grey$Module <- "Grey"
  grey$kME <- NA_real_
  grey$kME_p <- NA_real_
  grey$kME_q <- NA_real_
  module_members <- rbind(module_members, grey)
}
eigengene_df <- data.frame(
  GeneSampleID = sample_ids,
  eigengenes,
  TreatmentCode = meta$TreatmentCode,
  Year = meta$Year,
  Replicate = meta$Replicate,
  check.names = FALSE
)
# ---------------------------
# 7. Module eigengene analysis
# ---------------------------
association_rows <- list()
cell_factor <- traits$Cell
for (module in module_names) {
  for (trait_name in c("Yield", "Bacterial_state", "Metabolic_state")) {
    x <- eigengenes[, module]
    y <- traits[[trait_name]]
    overall <- cor.test(x, y, method = "spearman", exact = FALSE)
    adjusted <- robust_rank_association(x, y, cell_factor)
    association_rows[[length(association_rows) + 1]] <- data.frame(
      Module = module,
      Trait = trait_name,
      Overall_rho = unname(overall$estimate),
      Overall_p = overall$p.value,
      Adjusted_beta = adjusted["beta"],
      Adjusted_SE = adjusted["se"],
      Adjusted_CI_low = adjusted["ci_low"],
      Adjusted_CI_high = adjusted["ci_high"],
      Adjusted_p = adjusted["p"]
    )
  }
}
module_trait_assoc <- rbindlist(association_rows)
module_trait_assoc[, Overall_q := safe_bh(Overall_p)]
module_trait_assoc[, Adjusted_q := safe_bh(Adjusted_p)]
# ---------------------------
# 8. Module-by-trait association analysis
# ---------------------------
interaction_results <- list()
trajectory_results <- list()
for (module in module_names) {
  dat <- data.frame(
    ME = eigengenes[, module],
    TreatmentCode = factor(meta$TreatmentCode, levels = c("G", "H")),
    YearFactor = factor(meta$Year),
    Year = meta$Year
  )
  fit <- lm(ME ~ TreatmentCode * YearFactor, data = dat)
  wald <- robust_interaction_wald(fit, pattern = ":")
  interaction_results[[module]] <- data.frame(
    Module = module,
    Interaction_Wald = wald["statistic"],
    Interaction_df = wald["df"],
    Interaction_p = wald["p"]
  )
  for (year in sort(unique(dat$Year))) {
    sub <- dat[dat$Year == year, ]
    g <- sub$ME[sub$TreatmentCode == "G"]
    h <- sub$ME[sub$TreatmentCode == "H"]
    test <- t.test(h, g, var.equal = FALSE)
    trajectory_results[[length(trajectory_results) + 1]] <- data.frame(
      Module = module,
      Year = year,
      G_mean = mean(g),
      G_SE = sd(g) / sqrt(length(g)),
      H_mean = mean(h),
      H_SE = sd(h) / sqrt(length(h)),
      H_minus_G = mean(h) - mean(g),
      P = test$p.value
    )
  }
}
interaction_table <- rbindlist(interaction_results)
interaction_table[, Interaction_q := p.adjust(Interaction_p, method = "BH")]
module_trajectory <- rbindlist(trajectory_results)
module_trajectory[, Q_within_module := p.adjust(P, method = "BH"), by = Module]
# ---------------------------
# 9.
# ---------------------------
hub_functions <- module_members[
  Module != "Grey" & abs(kME) >= 0.80 & kME_q < 0.05
][order(Module, -abs(kME))]
# ---------------------------
# 10.
# ---------------------------
preservation_rows <- list()
for (left_year in sort(unique(meta$Year))) {
  keep <- meta$Year != left_year
  sub_cor <- cor(X_network[keep, , drop = FALSE], use = "pairwise.complete.obs")
  for (module in module_names) {
    idx <- which(module_assignment == module)
    within <- sub_cor[idx, idx, drop = FALSE]
    observed <- mean((1 + within[upper.tri(within)]) / 2, na.rm = TRUE)
    null_values <- replicate(n_preservation_permutations, {
      random_idx <- sample(seq_len(ncol(X_network)), length(idx), replace = FALSE)
      random_mat <- sub_cor[random_idx, random_idx, drop = FALSE]
      mean((1 + random_mat[upper.tri(random_mat)]) / 2, na.rm = TRUE)
    })
    z_pres <- (observed - mean(null_values)) / sd(null_values)
    p_perm <- (sum(null_values >= observed) + 1) / (length(null_values) + 1)
    preservation_rows[[length(preservation_rows) + 1]] <- data.frame(
      Left_out_year = left_year,
      Module = module,
      Mean_signed_similarity = observed,
      Z_preservation = z_pres,
      Permutation_p = p_perm
    )
  }
}
module_preservation <- rbindlist(preservation_rows)
# ---------------------------
# 11. M1
# ---------------------------
key_module <- "M1"
m1_keys <- names(module_assignment)[module_assignment == key_module]
common_m1_keys <- intersect(m1_keys, rownames(mg_all))
mg_m1 <- t(mg_all[common_m1_keys, , drop = FALSE])
mg_m1 <- apply(mg_m1, 2, z_score)
mg_m1 <- as.matrix(mg_m1)
rownames(mg_m1) <- sample_ids
m1_weights <- loading_list[[key_module]][common_m1_keys]
m1_weights <- m1_weights / sqrt(sum(m1_weights^2))
mg_m1_score <- as.numeric(mg_m1 %*% m1_weights)
mt_m1_score <- eigengenes[, key_module]
mg_mt_test <- cor.test(mg_m1_score, mt_m1_score, method = "pearson")
mg_adjusted <- robust_rank_association(mg_m1_score, mt_m1_score, cell_factor)
mg_validation <- data.frame(
  GeneSampleID = sample_ids,
  M1_MT_eigengene = mt_m1_score,
  M1_MG_projected_score = mg_m1_score,
  TreatmentCode = meta$TreatmentCode,
  Year = meta$Year
)
# ---------------------------
# 12. PLS
# SEMRT,
# ,
# ---------------------------
# SEMRT.
# ,.
fit_bacterial_state_fold <- function(train_idx, test_idx) {
  train_taxa <- taxa_mat[, train_idx, drop = FALSE]
  test_taxa  <- taxa_mat[, test_idx, drop = FALSE]
  positive <- train_taxa[is.finite(train_taxa) & train_taxa > 0]
  if (length(positive) == 0) stop("4,RT.")
  pseudo <- min(positive, na.rm = TRUE) / 2
  calc_state <- function(x) {
    lx <- log(x + pseudo)
    colMeans(lx[high_taxa, , drop = FALSE]) -
      colMeans(lx[low_taxa, , drop = FALSE])
  }
  list(
    train = as.numeric(calc_state(train_taxa)),
    test = as.numeric(calc_state(test_taxa)),
    pseudocount = pseudo
  )
}
# SEM
# log1p -> SD -> -> PCA PC1 -> .
fit_metabolic_state_fold <- function(train_idx, test_idx) {
  train_log <- met_log_sample[train_idx, all_metabolites, drop = FALSE]
  test_log  <- met_log_sample[test_idx, all_metabolites, drop = FALSE]
  mu <- colMeans(train_log, na.rm = TRUE)
  sigma <- apply(train_log, 2, function(x) {
    sqrt(mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
  })
  sigma[!is.finite(sigma) | sigma < 1e-10] <- 1
  train_z <- sweep(sweep(train_log, 2, mu, "-"), 2, sigma, "/")
  test_z  <- sweep(sweep(test_log,  2, mu, "-"), 2, sigma, "/")
  train_z[, negative_metabolites] <- -train_z[, negative_metabolites, drop = FALSE]
  test_z[, negative_metabolites]  <- -test_z[, negative_metabolites, drop = FALSE]
  pca <- prcomp(train_z, center = FALSE, scale. = FALSE)
  loading <- pca$rotation[, 1]
  train_score <- as.numeric(train_z %*% loading)
  test_score  <- as.numeric(test_z %*% loading)
  if (cor(train_score, rowMeans(train_z), use = "complete.obs") < 0) {
    loading <- -loading
    train_score <- -train_score
    test_score <- -test_score
  }
  list(
    train = train_score,
    test = test_score,
    loading = loading,
    variance = pca$sdev[1]^2 / sum(pca$sdev^2)
  )
}
# + RT + PLS.
# PC1.
fit_coordinated_state_fold <- function(train_idx, test_idx) {
  bac <- fit_bacterial_state_fold(train_idx, test_idx)
  met <- fit_metabolic_state_fold(train_idx, test_idx)
  train_traits <- cbind(
    Yield = as.numeric(yield_value[train_idx]),
    Bacterial_state = bac$train,
    Metabolic_state = met$train
  )
  test_traits <- cbind(
    Yield = as.numeric(yield_value[test_idx]),
    Bacterial_state = bac$test,
    Metabolic_state = met$test
  )
  mu <- colMeans(train_traits, na.rm = TRUE)
  sigma <- apply(train_traits, 2, sd, na.rm = TRUE)
  sigma[!is.finite(sigma) | sigma < 1e-10] <- 1
  train_scaled <- sweep(sweep(train_traits, 2, mu, "-"), 2, sigma, "/")
  test_scaled  <- sweep(sweep(test_traits,  2, mu, "-"), 2, sigma, "/")
  pca <- prcomp(train_scaled, center = FALSE, scale. = FALSE)
  loading <- pca$rotation[, 1]
  train_score <- as.numeric(train_scaled %*% loading)
  test_score <- as.numeric(test_scaled %*% loading)
  if (cor(train_score, rowMeans(train_scaled), use = "complete.obs") < 0) {
    loading <- -loading
    train_score <- -train_score
    test_score <- -test_score
  }
  list(
    train_score = train_score,
    test_score = test_score,
    loading = loading,
    variance = pca$sdev[1]^2 / sum(pca$sdev^2),
    metabolic_variance = met$variance,
    bacterial_pseudocount = bac$pseudocount
  )
}
fit_pls_model <- function(x_train, y_train, x_test, ncomp) {
  scaled <- safe_scale_train_test(x_train, x_test)
  x_names <- paste0("X", seq_len(ncol(scaled$train)))
  train_df <- as.data.frame(scaled$train)
  test_df  <- as.data.frame(scaled$test)
  names(train_df) <- x_names
  names(test_df) <- x_names
  train_df$Y <- y_train
  model <- pls::plsr(
    Y ~ .,
    data = train_df,
    ncomp = ncomp,
    scale = FALSE,
    validation = "none",
    method = "simpls"
  )
  prediction <- as.numeric(predict(model, newdata = test_df, ncomp = ncomp))
  list(model = model, prediction = prediction, scaled = scaled)
}
calculate_vip <- function(model, ncomp) {
 # X
  Tscore_raw <- pls::scores(model)
  if (is.null(Tscore_raw)) {
    stop("PLSscores,VIP.")
  }
  Tscore <- as.matrix(Tscore_raw)
 # PLSX
 # oscorespls/kernelplsloading.weights
 # simplsloading.weights,Xscoresprojection.
  W_raw <- model$loading.weights
  weight_source <- "loading.weights"
  if (is.null(W_raw)) {
    W_raw <- model$projection
    weight_source <- "projection"
  }
  if (is.null(W_raw)) {
    stop(
      "PLSloading.weights,projection,",
      "VIP.:",
      paste(names(model), collapse = ", ")
    )
  }
  W <- as.matrix(W_raw)
 # Y
  Yloadings_raw <- model$Yloadings
  if (is.null(Yloadings_raw)) {
    stop(
      "PLSYloadings,VIP.:",
      paste(names(model), collapse = ", ")
    )
  }
  Yloadings <- as.matrix(Yloadings_raw)
  available_components <- min(
    as.integer(ncomp),
    ncol(Tscore),
    ncol(W),
    max(dim(Yloadings))
  )
  if (
    length(available_components) != 1 ||
    !is.finite(available_components) ||
    available_components < 1
  ) {
    stop(
      "VIP:ncomp=", ncomp,
      ";scores=", paste(dim(Tscore), collapse = "×"),
      ";weights=", paste(dim(W), collapse = "×"),
      ";Yloadings=", paste(dim(Yloadings), collapse = "×")
    )
  }
  components <- seq_len(available_components)
  Tscore <- Tscore[, components, drop = FALSE]
  W <- W[, components, drop = FALSE]
 # plsYloadings1 × ncomp
  if (nrow(Yloadings) == 1 && ncol(Yloadings) >= available_components) {
    Q <- as.numeric(Yloadings[1, components, drop = TRUE])
  } else if (ncol(Yloadings) == 1 && nrow(Yloadings) >= available_components) {
    Q <- as.numeric(Yloadings[components, 1, drop = TRUE])
  } else if (ncol(Yloadings) >= available_components) {
    Q <- as.numeric(Yloadings[1, components, drop = TRUE])
  } else {
    stop(
      "YloadingsPLS:",
      paste(dim(Yloadings), collapse = " × "),
      ";=", available_components
    )
  }
  if (length(Q) != available_components || any(!is.finite(Q))) {
    stop("Y,VIP.")
  }
  ss_y <- colSums(Tscore^2) * Q^2
  total_ss_y <- sum(ss_y)
  if (!is.finite(total_ss_y) || total_ss_y <= 0) {
    stop(
      "PLSY,VIP:",
      paste(signif(ss_y, 5), collapse = ", ")
    )
  }
 # XVIP
  weight_norm <- colSums(W^2)
  if (any(!is.finite(weight_norm)) || any(weight_norm <= 0)) {
    stop(
      "PLS", weight_source,
      ",VIP."
    )
  }
  normalized_w2 <- sweep(W^2, 2, weight_norm, "/")
  p <- nrow(W)
  vip <- sqrt(
    p * rowSums(sweep(normalized_w2, 2, ss_y, "*")) / total_ss_y
  )
  feature_names <- rownames(W)
  if (is.null(feature_names) || length(feature_names) != length(vip)) {
    feature_names <- paste0("X", seq_along(vip))
  }
  names(vip) <- feature_names
  if (any(!is.finite(vip))) {
    stop("VIPNA,NaNInf.")
  }
  attr(vip, "weight_source") <- weight_source
  vip
}
nested_leave_one_year_pls <- function(X, years, model_name, max_comp = 4, collect_importance = FALSE) {
  unique_years <- sort(unique(years))
  prediction_rows <- list()
  importance_rows <- list()
  component_rows <- list()
  for (outer_year in unique_years) {
    outer_train <- years != outer_year
    outer_test <- years == outer_year
    state <- fit_coordinated_state_fold(
      outer_train,
      outer_test
    )
    y_train <- state$train_score
    y_test <- state$test_score
    train_years <- years[outer_train]
    max_allowed <- min(max_comp, ncol(X), sum(outer_train) - 2)
    candidate_components <- seq_len(max_allowed)
    inner_mse <- sapply(candidate_components, function(ncomp) {
      inner_pred <- rep(NA_real_, sum(outer_train))
      for (inner_year in sort(unique(train_years))) {
        inner_train <- train_years != inner_year
        inner_test <- train_years == inner_year
        fit <- fit_pls_model(
          X[outer_train, , drop = FALSE][inner_train, , drop = FALSE],
          y_train[inner_train],
          X[outer_train, , drop = FALSE][inner_test, , drop = FALSE],
          ncomp
        )
        inner_pred[inner_test] <- fit$prediction
      }
      mean((inner_pred - y_train)^2)
    })
    best_ncomp <- candidate_components[which.min(inner_mse)]
    outer_fit <- fit_pls_model(
      X[outer_train, , drop = FALSE],
      y_train,
      X[outer_test, , drop = FALSE],
      best_ncomp
    )
    prediction_rows[[length(prediction_rows) + 1]] <- data.frame(
      Model = model_name,
      Left_out_year = outer_year,
      GeneSampleID = sample_ids[outer_test],
      Observed_state = y_test,
      Predicted_state = outer_fit$prediction,
      Training_PC1_variance = state$variance,
      Training_metabolic_PC1_variance = state$metabolic_variance,
      Training_bacterial_pseudocount = state$bacterial_pseudocount,
      N_components = best_ncomp
    )
    component_rows[[length(component_rows) + 1]] <- data.frame(
      Model = model_name,
      Left_out_year = outer_year,
      N_components = best_ncomp,
      Inner_MSE = min(inner_mse),
      Training_PC1_variance = state$variance,
      Training_metabolic_PC1_variance = state$metabolic_variance,
      Training_bacterial_pseudocount = state$bacterial_pseudocount
    )
    if (collect_importance) {
      vip <- calculate_vip(outer_fit$model, best_ncomp)
      coefficient_array <- coef(
        outer_fit$model,
        ncomp = best_ncomp,
        intercept = FALSE
      )
      coefficient <- as.numeric(coefficient_array[, 1, 1])
      if (length(vip) != ncol(X)) {
        stop(
          "VIP:VIP=", length(vip),
          ";X=", ncol(X),
          ";=", outer_year
        )
      }
      if (length(coefficient) != ncol(X)) {
        stop(
          "PLS:=", length(coefficient),
          ";X=", ncol(X),
          ";=", outer_year
        )
      }
      importance_rows[[length(importance_rows) + 1]] <- data.frame(
        Model = model_name,
        Left_out_year = outer_year,
        FeatureKey = colnames(X),
        VIP = as.numeric(vip),
        Coefficient = coefficient
      )
    }
  }
  predictions <- rbindlist(prediction_rows)
  performance <- predictions[
    ,
    .(
      Pearson_r = cor(Observed_state, Predicted_state),
      CV_R2 = r2_score(Observed_state, Predicted_state),
      RMSE = sqrt(mean((Observed_state - Predicted_state)^2)),
      Median_training_PC1_variance = median(Training_PC1_variance),
      Median_components = median(N_components)
    ),
    by = Model
  ]
  list(
    predictions = predictions,
    performance = performance,
    components = rbindlist(component_rows),
    importance = if (collect_importance) rbindlist(importance_rows) else NULL
  )
}
X_mt <- t(mt_all)
X_mg <- t(mg_all)
X_modules <- eigengenes
X_mt <- X_mt[sample_ids, , drop = FALSE]
X_mg <- X_mg[sample_ids, , drop = FALSE]
X_modules <- X_modules[sample_ids, , drop = FALSE]
pls_mt <- nested_leave_one_year_pls(
  X_mt, meta$Year,
  model_name = "Metatranscriptome functional units",
  max_comp = max_pls_components,
  collect_importance = TRUE
)
pls_mg <- nested_leave_one_year_pls(
  X_mg, meta$Year,
  model_name = "Metagenome functional units",
  max_comp = max_pls_components,
  collect_importance = FALSE
)
pls_modules <- nested_leave_one_year_pls(
  X_modules, meta$Year,
  model_name = "Co-expression module eigengenes",
  max_comp = min(max_pls_components, ncol(X_modules)),
  collect_importance = FALSE
)
pls_predictions <- rbindlist(list(
  pls_mt$predictions, pls_mg$predictions, pls_modules$predictions
))
pls_performance <- rbindlist(list(
  pls_mt$performance, pls_mg$performance, pls_modules$performance
))
pls_components <- rbindlist(list(
  pls_mt$components, pls_mg$components, pls_modules$components
))
# ---------------------------
# 13. Integrated kME and VIP analysis
# ---------------------------
vip_summary <- pls_mt$importance[
  ,
  .(
    Mean_VIP = mean(VIP),
    VIP_frequency = mean(VIP > 1),
    Mean_coefficient = mean(Coefficient),
    Mean_abs_coefficient = mean(abs(Coefficient)),
    Sign_consistency = max(mean(Coefficient > 0), mean(Coefficient < 0))
  ),
  by = FeatureKey
]
m1_integrated <- merge(
  module_members[Module == "M1"],
  vip_summary,
  by = "FeatureKey",
  all.x = TRUE
)
m1_integrated[, Importance_score :=
                abs(kME) * Mean_VIP * Mean_abs_coefficient *
                pmax(VIP_frequency, 0.01) * pmax(Sign_consistency, 0.01)]
# M1
m1_category <- m1_integrated[
  is.finite(Importance_score),
  .(Weight = sum(Importance_score, na.rm = TRUE)),
  by = FunctionalModule
][order(-Weight)]
m1_category[, Contribution_percent := 100 * Weight / sum(Weight)]
# ---------------------------
# 14. :,,M1
# ---------------------------
trajectory_input <- data.frame(
  GeneSampleID = sample_ids,
  TreatmentCode = meta$TreatmentCode,
  Year = meta$Year,
  Yield = z_score(yield_value),
  Bacterial_state = z_score(bacterial_state),
  Metabolic_state = z_score(metabolic_state),
  M1_eigengene = z_score(eigengenes[, "M1"])
)
trajectory_long <- melt(
  as.data.table(trajectory_input),
  id.vars = c("GeneSampleID", "TreatmentCode", "Year"),
  variable.name = "Variable",
  value.name = "Value"
)
integrated_trajectory <- trajectory_long[
  ,
  .(
    G_mean = mean(Value[TreatmentCode == "G"]),
    H_mean = mean(Value[TreatmentCode == "H"]),
    Mowing_effect = mean(Value[TreatmentCode == "H"]) - mean(Value[TreatmentCode == "G"])
  ),
  by = .(Year, Variable)
]
trajectory_labels <- c(
  Yield = "Yield",
  Bacterial_state = "RT bacterial state",
  Metabolic_state = "Root metabolic state",
  M1_eigengene = "M1 eigengene"
)
# ---------------------------
# 15.
# ---------------------------
fwrite(metabolic_state_loadings, file.path(output_dir, "00_root_metabolic_state_PC1_loadings.csv"), bom = TRUE)
fwrite(power_table, file.path(output_dir, "01_soft_threshold_evaluation.csv"), bom = TRUE)
fwrite(module_summary, file.path(output_dir, "02_module_overview.csv"), bom = TRUE)
fwrite(data.table(
  FeatureKey = names(module_assignment),
  Module = module_assignment
), file.path(output_dir, "03_network_module_assignment.csv"), bom = TRUE)
fwrite(module_members, file.path(output_dir, "04_module_membership_kME.csv"), bom = TRUE)
fwrite(eigengene_df, file.path(output_dir, "05_module_eigengenes.csv"), bom = TRUE)
fwrite(module_trait_assoc, file.path(output_dir, "06_module_trait_associations.csv"), bom = TRUE)
fwrite(interaction_table, file.path(output_dir, "07_module_treatment_year_interaction.csv"), bom = TRUE)
fwrite(module_trajectory, file.path(output_dir, "08_module_year_trajectories.csv"), bom = TRUE)
fwrite(hub_functions, file.path(output_dir, "09_hub_functions.csv"), bom = TRUE)
fwrite(module_preservation, file.path(output_dir, "10_leave_one_year_module_preservation.csv"), bom = TRUE)
fwrite(mg_validation, file.path(output_dir, "11_M1_metagenome_validation_samples.csv"), bom = TRUE)
fwrite(data.table(
  M1_MT_features = length(m1_keys),
  M1_MG_matched_features = length(common_m1_keys),
  Overall_Pearson_r = unname(mg_mt_test$estimate),
  Overall_P = mg_mt_test$p.value,
  Adjusted_beta = mg_adjusted["beta"],
  Adjusted_P = mg_adjusted["p"],
  Adjusted_CI_low = mg_adjusted["ci_low"],
  Adjusted_CI_high = mg_adjusted["ci_high"]
), file.path(output_dir, "12_M1_metagenome_validation_summary.csv"), bom = TRUE)
fwrite(pls_predictions, file.path(output_dir, "13_PLS_leave_one_year_predictions.csv"), bom = TRUE)
fwrite(pls_performance, file.path(output_dir, "14_PLS_cross_year_performance.csv"), bom = TRUE)
fwrite(pls_components, file.path(output_dir, "15_PLS_component_selection.csv"), bom = TRUE)
fwrite(vip_summary, file.path(output_dir, "16_MT_function_unit_VIP_stability.csv"), bom = TRUE)
fwrite(m1_integrated, file.path(output_dir, "17_M1_kME_VIP_integration.csv"), bom = TRUE)
fwrite(m1_category, file.path(output_dir, "18_M1_functional_category_contribution.csv"), bom = TRUE)
fwrite(integrated_trajectory, file.path(output_dir, "19_integrated_temporal_trajectory.csv"), bom = TRUE)
# ---------------------------
# 16. : + PLS
# ---------------------------
# A Temporal trajectories shown as small multiples
trajectory_display_labels <- c(
  Yield = "Yield",
  Bacterial_state = "RT bacterial state",
  Metabolic_state = "Root metabolic state",
  M1_eigengene = "M1 eigengene"
)
trajectory_palette <- c(
  Yield = "#3C5488",
  Bacterial_state = "#00A087",
  Metabolic_state = "#E64B35",
  M1_eigengene = "#7E6148"
)
trajectory_plot <- copy(integrated_trajectory)
trajectory_plot$Display <- factor(
  trajectory_display_labels[trajectory_plot$Variable],
  levels = unname(trajectory_display_labels)
)
pA <- ggplot(
  trajectory_plot,
  aes(
    x = Year,
    y = Mowing_effect,
    group = Variable,
    colour = Variable
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    linetype = 2,
    colour = "grey62"
  ) +
  geom_line(linewidth = 0.85, lineend = "round") +
  geom_point(
    size = 2.4,
    shape = 21,
    fill = "white",
    stroke = 0.75
  ) +
  facet_wrap(~Display, ncol = 2) +
  scale_colour_manual(values = trajectory_palette) +
  scale_x_continuous(
    breaks = sort(unique(meta$Year)),
    expand = expansion(mult = c(0.03, 0.06))
  ) +
  labs(
    title = "Temporal alignment of mowing responses",
    subtitle = "Standardized mowing effects across yield, microbial, metabolic, and transcriptional states",
    x = "Stand age (years)",
    y = "Standardized mowing effect (H - G)"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7.5, hjust = 0),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.35, "cm"),
    plot.subtitle = element_text(size = 7.4)
  )
# B Cross-year prediction with identity and fitted relationships
pred_module <- copy(pls_modules$predictions)
perf_module <- pls_modules$performance[1]
pred_module$Left_out_year_numeric <- as.numeric(as.character(pred_module$Left_out_year))
prediction_limits <- range(
  c(pred_module$Observed_state, pred_module$Predicted_state),
  na.rm = TRUE
)
prediction_padding <- diff(prediction_limits) * 0.08
if (!is.finite(prediction_padding) || prediction_padding <= 0) {
  prediction_padding <- 0.2
}
prediction_limits <- prediction_limits + c(-prediction_padding, prediction_padding)
pB <- ggplot(
  pred_module,
  aes(
    x = Observed_state,
    y = Predicted_state,
    colour = Left_out_year_numeric
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    linewidth = 0.45,
    colour = "grey55"
  ) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.65,
    colour = "grey20",
    fill = "grey85",
    alpha = 0.45,
    inherit.aes = TRUE
  ) +
  geom_point(
    size = 2.5,
    alpha = 0.90,
    shape = 21,
    fill = "white",
    stroke = 0.85
  ) +
  scale_colour_viridis_c(
    option = "C",
    end = 0.90,
    breaks = sort(unique(pred_module$Left_out_year_numeric))
  ) +
  coord_equal(
    xlim = prediction_limits,
    ylim = prediction_limits,
    expand = FALSE
  ) +
  annotate(
    "label",
    x = prediction_limits[1],
    y = prediction_limits[2],
    hjust = 0,
    vjust = 1,
    size = 2.35,
    label.size = 0.25,
    label = sprintf(
      "r = %.2f\nCV R² = %.2f\nRMSE = %.2f\nPC1 variance = %.1f%%",
      perf_module$Pearson_r,
      perf_module$CV_R2,
      perf_module$RMSE,
      100 * perf_module$Median_training_PC1_variance
    )
  ) +
  labs(
    title = "Cross-year prediction from module eigengenes",
    subtitle = "Observed versus leave-one-year-out predicted coordinated states",
    x = "Observed coordinated state",
    y = "Predicted coordinated state",
    colour = "Left-out year"
  ) +
  guides(
    colour = guide_colourbar(
      title.position = "top",
      barheight = unit(2.1, "cm"),
      barwidth = unit(0.28, "cm")
    )
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(size = 7.4)
  )
# Shared functional-class palette for Panels C and E
functional_palette <- c(
  "K/P transporter validation" = "#E64B35",
  "Mineral-weathering candidate" = "#B09C00",
  "Phosphorus acquisition" = "#00A087",
  "Potassium acquisition and homeostasis" = "#4DBBD5",
  "Turnover-supporting functions" = "#DC6BE5",
  "Polyphenol-related functions" = "#7E6148",
  "CAZy resource acquisition and turnover" = "#8491B4"
)
# C Mini-bar + dot hybrid for top M1 features
m1_plot <- m1_integrated[
  is.finite(Mean_VIP) & is.finite(kME) & is.finite(Importance_score)
]
m1_plot[, abs_kME := abs(kME)]
m1_plot <- m1_plot[order(-Importance_score)]
top_n_display <- min(16, nrow(m1_plot))
m1_top <- copy(m1_plot[seq_len(top_n_display)])
m1_top$DisplayUnit <- stringr::str_wrap(m1_top$Unit, width = 24)
m1_top$DisplayUnit <- factor(m1_top$DisplayUnit, levels = rev(m1_top$DisplayUnit))
m1_top$kME_label <- sprintf("|kME| %.2f", m1_top$abs_kME)
x_max_c <- max(m1_top$Importance_score, na.rm = TRUE)
pC <- ggplot(
  m1_top,
  aes(x = Importance_score, y = DisplayUnit)
) +
  geom_col(
    aes(fill = FunctionalModule),
    width = 0.58,
    alpha = 0.28,
    show.legend = TRUE
  ) +
  geom_point(
    aes(size = Mean_VIP, fill = FunctionalModule),
    shape = 21,
    colour = "black",
    stroke = 0.30,
    alpha = 0.96
  ) +
  geom_text(
    aes(label = kME_label),
    hjust = 0,
    nudge_x = x_max_c * 0.035,
    size = 2.15,
    colour = "grey20",
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2.6, 6.0)) +
  scale_fill_manual(values = functional_palette, drop = FALSE) +
  coord_cartesian(
    xlim = c(0, x_max_c * 1.25),
    clip = "off"
  ) +
  labs(
    title = "Top M1 features ranked by integrated importance",
    subtitle = "Bar length indicates integrated importance; point size indicates mean VIP",
    x = "Integrated importance score",
    y = NULL,
    fill = "Functional class",
    size = "Mean VIP"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 6.4),
    legend.title = element_text(size = 7),
    axis.text.y = element_text(size = 6.7),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 58, 5.5, 5.5)
  )
# D Dumbbell comparison of predictive metrics
model_label_map <- c(
  "Metatranscriptome functional units" = "MT functional units",
  "Co-expression module eigengenes" = "Module eigengenes",
  "Metagenome functional units" = "MG functional units"
)
perf_plot <- copy(pls_performance)
perf_plot$Display <- unname(model_label_map[perf_plot$Model])
if (anyNA(perf_plot$Display)) {
  stop(
    "Unmatched PLS model labels: ",
    paste(perf_plot$Model[is.na(perf_plot$Display)], collapse = ", ")
  )
}
model_order <- c(
  "MT functional units",
  "Module eigengenes",
  "MG functional units"
)
perf_plot$Display <- factor(perf_plot$Display, levels = rev(model_order))
metric_min <- min(
  -0.05,
  perf_plot$CV_R2,
  perf_plot$Pearson_r,
  na.rm = TRUE
) - 0.05
metric_max <- max(
  1.00,
  perf_plot$CV_R2,
  perf_plot$Pearson_r,
  na.rm = TRUE
) + 0.08
pD <- ggplot(perf_plot, aes(y = Display)) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    linetype = 2,
    colour = "grey60"
  ) +
  geom_segment(
    aes(x = CV_R2, xend = Pearson_r, yend = Display),
    linewidth = 1.05,
    colour = "grey75",
    lineend = "round"
  ) +
  geom_point(
    aes(x = CV_R2),
    shape = 21,
    size = 3.4,
    fill = "#3C5488",
    colour = "black",
    stroke = 0.45
  ) +
  geom_point(
    aes(x = Pearson_r),
    shape = 22,
    size = 3.2,
    fill = "#E64B35",
    colour = "black",
    stroke = 0.45
  ) +
  geom_text(
    aes(x = CV_R2, label = sprintf("%.2f", CV_R2)),
    nudge_y = 0.17,
    size = 2.25,
    colour = "#3C5488"
  ) +
  geom_text(
    aes(x = Pearson_r, label = sprintf("%.2f", Pearson_r)),
    nudge_y = -0.17,
    size = 2.25,
    colour = "#E64B35"
  ) +
  geom_text(
    aes(
      x = metric_max - 0.01,
      label = sprintf("RMSE %.2f", RMSE)
    ),
    hjust = 1,
    size = 2.20,
    colour = "grey25"
  ) +
  scale_x_continuous(
    limits = c(metric_min, metric_max),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(
    title = "Cross-year predictive performance",
    subtitle = "Connected estimates show CV R² and Pearson r for each representation",
    x = "Predictive metric",
    y = NULL
  ) +
  annotate(
    "point",
    x = metric_min + 0.03 * (metric_max - metric_min),
    y = Inf,
    shape = 21,
    size = 2.8,
    fill = "#3C5488",
    colour = "black",
    stroke = 0.4
  ) +
  annotate(
    "text",
    x = metric_min + 0.07 * (metric_max - metric_min),
    y = Inf,
    label = "CV R²",
    hjust = 0,
    vjust = 1.4,
    size = 2.25
  ) +
  annotate(
    "point",
    x = metric_min + 0.28 * (metric_max - metric_min),
    y = Inf,
    shape = 22,
    size = 2.7,
    fill = "#E64B35",
    colour = "black",
    stroke = 0.4
  ) +
  annotate(
    "text",
    x = metric_min + 0.32 * (metric_max - metric_min),
    y = Inf,
    label = "Pearson r",
    hjust = 0,
    vjust = 1.4,
    size = 2.25
  ) +
  coord_cartesian(clip = "off") +
  theme(
    axis.text.y = element_text(size = 7.2),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 5.5, 5.5, 5.5),
    plot.subtitle = element_text(size = 7.2)
  )
# E Ranked contribution of functional classes
category_label_map <- c(
  "Turnover-supporting functions" = "Turnover",
  "Phosphorus acquisition" = "P acquisition",
  "Mineral-weathering candidate" = "Mineral weathering",
  "K/P transporter validation" = "K/P transport",
  "Potassium acquisition and homeostasis" = "K acquisition and homeostasis",
  "Polyphenol-related functions" = "Polyphenol-related",
  "CAZy resource acquisition and turnover" = "CAZy resource acquisition and turnover"
)
category_plot <- copy(m1_category)
category_plot$Display <- unname(
  category_label_map[category_plot$FunctionalModule]
)
category_plot$Display[is.na(category_plot$Display)] <-
  category_plot$FunctionalModule[is.na(category_plot$Display)]
category_plot <- category_plot[order(Contribution_percent)]
category_plot$Display <- factor(
  category_plot$Display,
  levels = category_plot$Display
)
e_max <- max(category_plot$Contribution_percent, na.rm = TRUE)
pE <- ggplot(
  category_plot,
  aes(
    x = Contribution_percent,
    y = Display,
    fill = FunctionalModule
  )
) +
  geom_col(width = 0.60, alpha = 0.82, show.legend = FALSE) +
  geom_point(
    aes(x = Contribution_percent),
    shape = 21,
    size = 2.8,
    fill = "white",
    colour = "black",
    stroke = 0.45,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", Contribution_percent)),
    hjust = 0,
    nudge_x = e_max * 0.035,
    size = 2.25,
    colour = "grey20"
  ) +
  scale_fill_manual(values = functional_palette, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, e_max * 1.20),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Functional composition of the stable M1 signal",
    subtitle = "Integrated contribution of functional classes",
    x = "Contribution (%)",
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 6.8),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(size = 7.2)
  )
pDE <- pD / pE + plot_layout(heights = c(0.95, 1.05))
main_figure <- (
  (pA | pB) /
    (pC | pDE)
) +
  plot_layout(
    widths = c(1.30, 1.00),
    heights = c(0.92, 1.28)
  ) +
  plot_annotation(
    tag_levels = "A",
    title = "Integrated analysis of RT functional co-expression and cross-year prediction",
    subtitle = "Temporal concordance, module-level prediction, and functional prioritization",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.8),
      plot.tag = element_text(face = "bold", size = 11)
    )
  )
ggsave(
  file.path(output_dir, "Figure4_RT_coexpression_PLS_integrated.png"),
  main_figure,
  width = 14.2,
  height = 10.8,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(output_dir, "Figure4_RT_coexpression_PLS_integrated.pdf"),
  main_figure,
  width = 14.2,
  height = 10.8
)
ggsave(
  file.path(output_dir, "Figure4_RT_coexpression_PLS_integrated.svg"),
  main_figure,
  width = 14.2,
  height = 10.8,
  device = svglite::svglite,
  bg = "white"
)
ggsave(
  file.path(output_dir, "Figure4_RT_coexpression_PLS_integrated_600dpi.tiff"),
  main_figure,
  width = 14.2,
  height = 10.8,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
# Export single panels for custom figure assembly
save_single_panel <- function(plot_obj, filename_base, width, height) {
  ggsave(
    file.path(single_panel_dir, paste0(filename_base, ".png")),
    plot_obj,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(single_panel_dir, paste0(filename_base, ".pdf")),
    plot_obj,
    width = width,
    height = height
  )
  ggsave(
    file.path(single_panel_dir, paste0(filename_base, ".svg")),
    plot_obj,
    width = width,
    height = height,
    device = svglite::svglite,
    bg = "white"
  )
  ggsave(
    file.path(single_panel_dir, paste0(filename_base, "_600dpi.tiff")),
    plot_obj,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}
save_single_panel(pA, "Panel_A_Temporal_trajectories", 7.6, 5.8)
save_single_panel(pB, "Panel_B_Cross_year_prediction", 6.3, 5.8)
save_single_panel(pC, "Panel_C_M1_feature_ranking", 7.8, 6.8)
save_single_panel(pD, "Panel_D_Predictive_performance", 6.3, 4.2)
save_single_panel(pE, "Panel_E_Functional_contribution", 6.3, 4.8)
# ---------------------------
# 17.
# ---------------------------
# S1 +
dend <- as.dendrogram(tree)
dend_data <- ggdendro::dendro_data(dend, type = "rectangle")
leaf_order <- order.dendrogram(dend)
strip_df <- data.frame(
  x = seq_along(leaf_order),
  y = -0.035,
  Module = module_assignment[tree$order],
  stringsAsFactors = FALSE
)
pS1 <- ggplot() +
  geom_segment(
    data = dend_data$segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.25,
    colour = "grey25"
  ) +
  geom_hline(yintercept = tree_cut_height, linetype = 2, linewidth = 0.4, colour = "grey55") +
  geom_tile(
    data = strip_df,
    aes(x = x, y = y, fill = Module),
    width = 1,
    height = 0.035
  ) +
  scale_fill_manual(values = module_colors) +
  coord_cartesian(ylim = c(-0.07, max(dend_data$segments$y) * 1.03), clip = "off") +
  labs(
    title = "Functional dendrogram and module assignment",
    subtitle = "Hierarchical clustering of functional units based on TOM dissimilarity",
    x = NULL,
    y = "1 - TOM",
    fill = "Module"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )
# Figure S2 diagnostic panels
heat <- melt(
  module_trait_assoc,
  id.vars = c("Module", "Trait"),
  measure.vars = c("Overall_rho", "Adjusted_beta"),
  variable.name = "Association_type",
  value.name = "Effect"
)
heat <- merge(
  heat,
  melt(
    module_trait_assoc,
    id.vars = c("Module", "Trait"),
    measure.vars = c("Overall_q", "Adjusted_q"),
    variable.name = "Q_type",
    value.name = "Q"
  )[, .(Module, Trait, Q_type, Q)],
  by = c("Module", "Trait"),
  allow.cartesian = TRUE
)
heat <- heat[
  (Association_type == "Overall_rho" & Q_type == "Overall_q") |
    (Association_type == "Adjusted_beta" & Q_type == "Adjusted_q")
]
heat$Association_type <- factor(
  heat$Association_type,
  levels = c("Overall_rho", "Adjusted_beta"),
  labels = c("Overall Spearman rho", "Adjusted beta")
)
heat$Label <- sprintf("%.2f%s", heat$Effect, ifelse(heat$Q < 0.05, "**", ifelse(heat$Q < 0.10, "†", "")))
pS2 <- ggplot(heat, aes(x = Trait, y = Module, fill = Effect)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = 2.6) +
  facet_wrap(~Association_type, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
  scale_x_discrete(labels = c(
    Yield = "Yield",
    Bacterial_state = "RT bacterial state",
    Metabolic_state = "Root metabolic state"
  )) +
  labs(
    title = "Module eigengene-trait associations",
    subtitle = "** q < 0.05; † q < 0.10",
    x = NULL,
    y = NULL,
    fill = "Effect"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
# S3 M1
m1_traj <- module_trajectory[Module == "M1"]
m1_long <- rbind(
  m1_traj[, .(Year, Treatment = "Control", Mean = G_mean, SE = G_SE)],
  m1_traj[, .(Year, Treatment = "Mowing", Mean = H_mean, SE = H_SE)]
)
pS3 <- ggplot(m1_long, aes(x = Year, y = Mean, colour = Treatment, shape = Treatment)) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.25, linewidth = 0.45) +
  scale_x_continuous(breaks = sort(unique(meta$Year))) +
  scale_colour_manual(values = c(Control = "#777777", Mowing = "#0072B2")) +
  labs(
    title = "Treatment-specific trajectory of the M1 eigengene",
    x = "Stand age (years)",
    y = "M1 eigengene",
    colour = NULL,
    shape = NULL
  ) +
  theme(legend.position = "top")
# S4
pS4 <- ggplot(
  module_preservation,
  aes(x = Left_out_year, y = Z_preservation, colour = Module, shape = Module)
) +
  geom_hline(yintercept = 10, linetype = 2, linewidth = 0.4, colour = "grey55") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = module_colors[module_names]) +
  scale_x_continuous(breaks = sort(unique(meta$Year))) +
  labs(
    title = "Leave-one-year-out module preservation",
    x = "Left-out year",
    y = "Preservation Z",
    colour = NULL,
    shape = NULL
  ) +
  theme(legend.position = "top")
# S5 M1
mg_validation_plot <- as.data.table(copy(mg_validation))
mg_validation_plot[, Treatment := factor(
  TreatmentCode,
  levels = c("G", "H"),
  labels = c("Control", "Mowing")
)]
mg_validation_plot[, Year_numeric := as.numeric(as.character(Year))]
mg_annotation <- sprintf(
  paste0(
    "Matched functions = %d/%d\n",
    "Overall Pearson r = %.3f, P = %.2e\n",
    "Adjusted beta = %.3f (95%% CI %.3f to %.3f), P = %.3f"
  ),
  length(common_m1_keys),
  length(m1_keys),
  unname(mg_mt_test$estimate),
  mg_mt_test$p.value,
  unname(mg_adjusted["beta"]),
  unname(mg_adjusted["ci_low"]),
  unname(mg_adjusted["ci_high"]),
  unname(mg_adjusted["p"])
)
mg_x_range <- range(
  mg_validation_plot$M1_MG_projected_score,
  na.rm = TRUE
)
mg_y_range <- range(
  mg_validation_plot$M1_MT_eigengene,
  na.rm = TRUE
)
pS5 <- ggplot(
  mg_validation_plot,
  aes(
    x = M1_MG_projected_score,
    y = M1_MT_eigengene
  )
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.70,
    colour = "grey20",
    fill = "grey82",
    alpha = 0.45
  ) +
  geom_point(
    aes(
      colour = Year_numeric,
      shape = Treatment
    ),
    size = 2.7,
    fill = "white",
    stroke = 0.85,
    alpha = 0.92
  ) +
  scale_colour_viridis_c(
    option = "C",
    end = 0.90,
    breaks = sort(unique(mg_validation_plot$Year_numeric))
  ) +
  scale_shape_manual(
    values = c(Control = 21, Mowing = 24)
  ) +
  annotate(
    "label",
    x = mg_x_range[1],
    y = mg_y_range[2],
    label = mg_annotation,
    hjust = 0,
    vjust = 1,
    size = 2.45,
    label.size = 0.25,
    fill = alpha("white", 0.90)
  ) +
  labs(
    title = "Metagenomic potential partially tracks the M1 transcriptional state",
    subtitle = "Overall association and treatment-by-year-adjusted robust rank regression",
    x = "M1 metagenomic functional-potential score",
    y = "M1 metatranscriptomic eigengene",
    colour = "Stand age",
    shape = NULL
  ) +
  guides(
    colour = guide_colourbar(
      title.position = "top",
      barheight = unit(2.0, "cm"),
      barwidth = unit(0.28, "cm")
    )
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    panel.grid.major = element_line(
      colour = "grey93",
      linewidth = 0.30
    ),
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(size = 7.4)
  )
supp_figure <- (
  (pS1 | pS2) /
    (pS3 | pS4) /
    pS5
) +
  plot_layout(heights = c(1.05, 1.00, 0.95)) +
  plot_annotation(tag_levels = "A")
ggsave(
  file.path(output_dir, "FigureS2_RT_coexpression_diagnostics.png"),
  supp_figure, width = 13.2, height = 13.2, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "FigureS2_RT_coexpression_diagnostics.pdf"),
  supp_figure, width = 13.2, height = 13.2
)
ggsave(
  file.path(output_dir, "FigureS2_RT_coexpression_diagnostics.svg"),
  supp_figure,
  width = 13.2,
  height = 13.2,
  device = svglite::svglite,
  bg = "white"
)
ggsave(
  file.path(output_dir, "FigureS2_RT_coexpression_diagnostics_600dpi.tiff"),
  supp_figure,
  width = 13.2,
  height = 13.2,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
# Export the five supplementary panels separately
save_supplementary_panel <- function(plot_obj, filename_base, width, height) {
  ggsave(
    file.path(supplementary_panel_dir, paste0(filename_base, ".png")),
    plot_obj,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(supplementary_panel_dir, paste0(filename_base, ".pdf")),
    plot_obj,
    width = width,
    height = height
  )
  ggsave(
    file.path(supplementary_panel_dir, paste0(filename_base, ".svg")),
    plot_obj,
    width = width,
    height = height,
    device = svglite::svglite,
    bg = "white"
  )
  ggsave(
    file.path(supplementary_panel_dir, paste0(filename_base, "_600dpi.tiff")),
    plot_obj,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}
save_supplementary_panel(
  pS1,
  "FigS2A_Functional_dendrogram",
  width = 8.2,
  height = 5.0
)
save_supplementary_panel(
  pS2,
  "FigS2B_Module_trait_associations",
  width = 8.2,
  height = 5.2
)
save_supplementary_panel(
  pS3,
  "FigS2C_M1_treatment_trajectory",
  width = 6.4,
  height = 4.8
)
save_supplementary_panel(
  pS4,
  "FigS2D_Module_preservation",
  width = 6.4,
  height = 4.8
)
save_supplementary_panel(
  pS5,
  "FigS2E_M1_MG_MT_validation",
  width = 7.2,
  height = 5.2
)
# ---------------------------
# ---------------------------
# Final reproducibility record
# ---------------------------
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo_Fig4_FigS2.txt")
)
message("Figure 4 and Figure S2 analyses completed.")
