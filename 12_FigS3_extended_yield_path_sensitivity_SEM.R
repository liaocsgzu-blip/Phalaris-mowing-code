# ======================================================================
# SEM from raw files + two model comparison + shared-variance analysis
# :M1,RT,,8,
# ,.
# Expected raw files in the working directory
# sample_mapping*.csv
# metatranscriptome_selected_functions*.csv
# rt_species_abundance*.csv
# plant_traits_yield*.csv
# root_metabolome_summary*.csv
# carbon_use_rs_bs*.csv
# Model A: no RT bacterial state -> Root activity branch
# Model B: add RT bacterial state -> Root activity branch
# BOTH models additionally include direct paths to Yield from
# Mowing, RT bacterial state, Microbial turnover, M1 functional module,
# Root metabolic module, and Root activity.
# ,RT,,M1,.
# Main shared-variance plot
# RT bacterial state + Root metabolic module -> Root activity
# after partialling out Year + Year^2
# ======================================================================
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260730)
# ----------------------------------------------------------------------
# 0. Packages
# ----------------------------------------------------------------------
pkgs <- c(
  "readr", "dplyr", "tidyr", "stringr", "tibble",
  "ggplot2", "sandwich", "lmtest", "svglite"
)
missing_pkgs <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Please install required packages first:\n",
    paste(missing_pkgs, collapse = ", "),
    "\n\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))"
  )
}
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(sandwich)
library(lmtest)
# ----------------------------------------------------------------------
# 1. Paths and helper functions
# ----------------------------------------------------------------------
ROOT <- getwd()
OUTDIR <- file.path(ROOT, "FigS3_extended_SEM_output")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
find_file <- function(pattern) {
  x <- list.files(
    ROOT,
    pattern = pattern,
    full.names = TRUE
  )
  if (length(x) == 0) {
    stop("Cannot find file matching: ", pattern)
  }
  x[1]
}
f_map <- find_file("^sample_mapping.*\\.csv$")
f_mt  <- find_file("^metatranscriptome_selected_functions.*\\.csv$")
f_16s <- find_file("^rt_species_abundance.*\\.csv$")
f_plant <- find_file("^plant_traits_yield.*\\.csv$")
f_metab <- find_file("^root_metabolome_summary.*\\.csv$")
f_carbon <- find_file("^carbon_use_rs_bs.*\\.csv$")
cat("Files detected:\n")
cat("Map: ", f_map, "\n")
cat("Metatranscriptome: ", f_mt, "\n")
cat("16S: ", f_16s, "\n")
cat("Plant: ", f_plant, "\n")
cat("Metabolome: ", f_metab, "\n")
cat("Carbon use: ", f_carbon, "\n\n")
meta <- read_csv(f_map, show_col_types = FALSE)
mt <- read_csv(f_mt, show_col_types = FALSE)
sp <- read_csv(f_16s, show_col_types = FALSE)
plant <- read_csv(f_plant, show_col_types = FALSE)
metab <- read_csv(f_metab, show_col_types = FALSE)
carbon <- read_csv(f_carbon, show_col_types = FALSE)
rt_samples <- meta$GeneSampleID
# Population-SD standardization, matching the analysis used here.
zpop <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  s <- sqrt(mean((x - mu)^2, na.rm = TRUE))
  (x - mu) / s
}
hc3_table <- function(model) {
  ct <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  tibble(
    Term = rownames(ct),
    Estimate = ct[, 1],
    SE_HC3 = ct[, 2],
    t = ct[, 3],
    P = ct[, 4]
  )
}
get_hc3 <- function(model, term) {
  x <- hc3_table(model)
  x[x$Term == term, , drop = FALSE]
}
# ----------------------------------------------------------------------
# 2. Rebuild M1 from the raw RT metatranscriptome
# Exact co-expression workflow used previously
# prevalence >= 30%; database-wise relative abundance; database-
# specific pseudocount; CLR; top 250 by MAD; signed adjacency power 10
# TOM; average-linkage clustering; cut height 0.905; minimum module 12.
# ----------------------------------------------------------------------
PREVALENCE <- 0.30
TOP_N <- 250L
POWER <- 10
CUT_HEIGHT <- 0.905
MIN_MODULE_SIZE <- 12L
missing_samples <- setdiff(rt_samples, names(mt))
if (length(missing_samples) > 0) {
  stop(
    "Metatranscriptome file is missing RT samples:\n",
    paste(missing_samples, collapse = ", ")
  )
}
mt[rt_samples] <- lapply(
  mt[rt_samples],
  function(x) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- 0
    x
  }
)
extract_pipe_field <- function(x, i) {
  z <- str_split(x, "\\s*\\|\\s*")
  vapply(
    z,
    function(parts) {
      if (length(parts) >= i) str_trim(parts[i]) else ""
    },
    character(1)
  )
}
make_unit <- function(Database, FeatureID, Annotation) {
  if (Database == "KEGG") {
    gene <- str_match(Annotation, "\\|\\s*([^;|]+)")[, 2]
    gene <- str_trim(gene)
    if (is.na(gene) || gene == "") return(as.character(FeatureID))
    return(paste0(FeatureID, "|", gene))
  }
  if (Database == "CAZy") {
    return(as.character(FeatureID))
  }
  if (Database == "Pcyc") {
    p3 <- extract_pipe_field(Annotation, 3)
    p4 <- extract_pipe_field(Annotation, 4)
    p5 <- extract_pipe_field(Annotation, 5)
    return(paste(p3, p4, p5, sep = "|"))
  }
  if (Database == "Polyphenol") {
    p3 <- extract_pipe_field(Annotation, 3)
    p4 <- extract_pipe_field(Annotation, 4)
    p5 <- extract_pipe_field(Annotation, 5)
    return(paste(p3, p4, p5, sep = "|"))
  }
  if (Database == "TCDB") {
    tcid <- str_extract(
      Annotation,
      "\\b[0-9]+\\.[A-Z]\\.[0-9]+(?:\\.[0-9]+)+\\b"
    )
    gene <- str_match(
      Annotation,
      "\\bGN=([A-Za-z0-9_.-]+)\\b"
    )[, 2]
    if (!is.na(gene) && gene != "" && !is.na(tcid) && tcid != "") {
      return(paste0(gene, "|", tcid))
    }
    if (!is.na(tcid) && tcid != "") return(tcid)
    if (!is.na(gene) && gene != "") return(gene)
    return(Annotation)
  }
  as.character(FeatureID)
}
mt$Unit <- mapply(
  make_unit,
  mt$Database,
  mt$FeatureID,
  mt$Annotation,
  USE.NAMES = FALSE
)
mt_agg <- mt %>%
  select(
    Database,
    FunctionalModule = Module,
    Unit,
    all_of(rt_samples)
  ) %>%
  group_by(Database, FunctionalModule, Unit) %>%
  summarise(
    across(all_of(rt_samples), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )
abund_mat <- as.matrix(mt_agg[, rt_samples])
prevalence <- rowMeans(abund_mat > 0, na.rm = TRUE)
mt_keep <- mt_agg[
  prevalence >= PREVALENCE,
  ,
  drop = FALSE
]
cat("M1 rebuild: retained after prevalence filter = ", nrow(mt_keep), "\n")
if (nrow(mt_keep) != 385) {
  warning(
    "Previous analysis retained 385 functional units after prevalence filtering; ",
    "current retained count = ", nrow(mt_keep),
    ". Check raw-file version if unexpected."
  )
}
clr_blocks <- list()
pseudo_table <- list()
for (db in unique(mt_keep$Database)) {
  block <- mt_keep %>% filter(Database == db)
  X <- as.matrix(block[, rt_samples])
  sample_totals <- colSums(X, na.rm = TRUE)
  sample_totals[
    !is.finite(sample_totals) | sample_totals <= 0
  ] <- 1
  rel <- sweep(X, 2, sample_totals, "/")
  positive_values <- rel[is.finite(rel) & rel > 0]
  pseudo <- 0.5 * unname(
    quantile(
      positive_values,
      probs = 0.001,
      na.rm = TRUE,
      type = 7
    )
  )
  pseudo_table[[db]] <- tibble(
    Database = db,
    Pseudocount = pseudo
  )
  log_rel <- log(rel + pseudo)
  clr <- sweep(
    log_rel,
    2,
    colMeans(log_rel, na.rm = TRUE),
    "-"
  )
  rownames(clr) <- paste(
    block$Database,
    block$FunctionalModule,
    block$Unit,
    sep = "||"
  )
  clr_blocks[[db]] <- clr
}
clr_feature_sample <- do.call(rbind, clr_blocks)
sds <- apply(clr_feature_sample, 1, sd, na.rm = TRUE)
keep_var <- is.finite(sds) & sds > 1e-12
clr_feature_sample <- clr_feature_sample[keep_var, , drop = FALSE]
mads <- apply(
  clr_feature_sample,
  1,
  function(v) {
    v <- as.numeric(v)
    stats::mad(
      v,
      center = stats::median(v, na.rm = TRUE),
      constant = 1,
      na.rm = TRUE
    )
  }
)
ord <- order(mads, decreasing = TRUE)
top_ids <- rownames(clr_feature_sample)[
  ord[seq_len(min(TOP_N, length(ord)))]
]
net_feature_sample <- clr_feature_sample[top_ids, , drop = FALSE]
# samples x features
X <- t(net_feature_sample)
# Feature-wise z standardization
Xz <- scale(X, center = TRUE, scale = TRUE)
if (any(!is.finite(Xz))) {
  stop("Non-finite values after feature standardization.")
}
R <- cor(
  Xz,
  use = "pairwise.complete.obs",
  method = "pearson"
)
adj <- ((1 + R) / 2)^POWER
diag(adj) <- 0
connectivity <- rowSums(adj, na.rm = TRUE)
L <- adj %*% adj
denominator <- outer(connectivity, connectivity, pmin) + 1 - adj
TOM <- (L + adj) / denominator
diag(TOM) <- 1
TOM[!is.finite(TOM)] <- 0
dissTOM <- 1 - TOM
diag(dissTOM) <- 0
tree <- hclust(
  as.dist(dissTOM),
  method = "average"
)
raw_labels <- cutree(tree, h = CUT_HEIGHT)
raw_sizes <- table(raw_labels)
small_clusters <- as.integer(
  names(raw_sizes[raw_sizes < MIN_MODULE_SIZE])
)
clean_labels <- raw_labels
clean_labels[clean_labels %in% small_clusters] <- 0
formal_sizes <- sort(
  table(clean_labels[clean_labels != 0]),
  decreasing = TRUE
)
label_to_module <- setNames(
  paste0("M", seq_along(formal_sizes)),
  names(formal_sizes)
)
Module <- ifelse(
  clean_labels == 0,
  "Grey",
  unname(label_to_module[as.character(clean_labels)])
)
module_table <- tibble(
  FeatureKey = colnames(Xz),
  RawCluster = raw_labels,
  CleanCluster = clean_labels,
  Module = Module,
  MAD = mads[colnames(Xz)],
  Connectivity = connectivity
)
module_sizes <- module_table %>%
  count(Module, name = "N") %>%
  arrange(factor(Module, levels = c("M1","M2","M3","M4","Grey")))
write_csv(
  module_sizes,
  file.path(OUTDIR, "01_M1_module_sizes.csv")
)
cat("\nModule sizes:\n")
print(module_sizes)
expected <- c(M1 = 93, M2 = 73, M3 = 15, M4 = 12, Grey = 57)
observed <- setNames(module_sizes$N, module_sizes$Module)
if (
  !all(
    observed[intersect(names(expected), names(observed))] ==
      expected[intersect(names(expected), names(observed))]
  )
) {
  warning(
    "Module sizes differ from expected M1=93, M2=73, M3=15, M4=12, Grey=57."
  )
}
formal_modules <- module_sizes$Module[module_sizes$Module != "Grey"]
ME_list <- list()
for (mod in formal_modules) {
  ids <- module_table$FeatureKey[module_table$Module == mod]
  Xm <- Xz[, ids, drop = FALSE]
  pc <- prcomp(
    Xm,
    center = FALSE,
    scale. = FALSE
  )
  score <- pc$x[, 1]
 # Orient PC1 toward higher average standardized module expression.
  if (cor(score, rowMeans(Xm), use = "complete.obs") < 0) {
    score <- -score
    pc$rotation[, 1] <- -pc$rotation[, 1]
  }
  pc1_var <- summary(pc)$importance[2, 1]
  ME_list[[mod]] <- tibble(
    GeneSampleID = rownames(Xm),
    Module = mod,
    Eigengene = score,
    PC1_variance_explained = pc1_var
  )
}
ME_long <- bind_rows(ME_list)
M1 <- ME_long %>%
  filter(Module == "M1") %>%
  select(
    GeneSampleID,
    M1 = Eigengene,
    M1_PC1_variance = PC1_variance_explained
  ) %>%
  mutate(GeneSampleID = as.character(GeneSampleID)) %>%
  arrange(match(GeneSampleID, rt_samples))
if (nrow(M1) != length(rt_samples) ||
    anyNA(match(rt_samples, M1$GeneSampleID))) {
  stop(
    "M1 sample IDs do not match metadata / M1ID.\n",
    "Missing from M1 / M1: ",
    paste(setdiff(rt_samples, M1$GeneSampleID), collapse = ", ")
  )
}
cat(
  "\nM1 size = ",
  sum(module_table$Module == "M1"),
  "; M1 PC1 variance = ",
  round(100 * unique(M1$M1_PC1_variance), 1),
  "%\n",
  sep = ""
)
write_csv(
  M1,
  file.path(OUTDIR, "02_M1_module_eigengene_39samples.csv")
)
# ----------------------------------------------------------------------
# 3. Build RT bacterial state from the four pre-selected taxa
# ----------------------------------------------------------------------
high_taxa <- c(
  "Staphylococcus capitis",
  "Pseudomonas syncyanea"
)
low_taxa <- c(
  "Luteitalea pratensis",
  "Paraflavitalea soli"
)
four_taxa <- c(high_taxa, low_taxa)
sp4 <- sp %>%
  filter(Species %in% four_taxa)
if (nrow(sp4) != 4) {
  stop(
    "Cannot find all four key taxa. Found: ",
    paste(sp4$Species, collapse = ", ")
  )
}
sp_mat <- as.matrix(sp4[, rt_samples])
storage.mode(sp_mat) <- "numeric"
positive_sp <- sp_mat[is.finite(sp_mat) & sp_mat > 0]
# Global pseudocount used in the previous analysis
# half of the smallest positive abundance across all four taxa/samples.
sp_pseudo <- 0.5 * min(positive_sp)
log_sp <- log(sp_mat + sp_pseudo)
rownames(log_sp) <- sp4$Species
bacterial_balance <- (
  colMeans(log_sp[high_taxa, , drop = FALSE]) -
  colMeans(log_sp[low_taxa, , drop = FALSE])
)
bacteria_df <- tibble(
  GeneSampleID = as.character(colnames(log_sp)),
  BacterialBalance = as.numeric(bacterial_balance)
) %>%
  arrange(match(GeneSampleID, rt_samples))
if (nrow(bacteria_df) != length(rt_samples) ||
    anyNA(match(rt_samples, bacteria_df$GeneSampleID))) {
  stop(
    "RT bacterial-state sample IDs do not match metadata / RTID."
  )
}
# ----------------------------------------------------------------------
# 4. Build microbial turnover state from RS + BS turnover PC1
# ----------------------------------------------------------------------
turnover_wide <- carbon %>%
  select(
    GeneSampleID,
    SoilInterface,
    Turnover = `Turnover rate`
  ) %>%
  pivot_wider(
    names_from = SoilInterface,
    values_from = Turnover,
    names_prefix = "Turnover_"
  )
turnover_wide <- turnover_wide %>%
  mutate(GeneSampleID = as.character(GeneSampleID)) %>%
  arrange(match(GeneSampleID, rt_samples))
missing_turnover_ids <- setdiff(rt_samples, turnover_wide$GeneSampleID)
if (length(missing_turnover_ids) > 0) {
  stop(
    "Turnover table is missing RT samples:\n",
    paste(missing_turnover_ids, collapse = ", ")
  )
}
turnover_wide <- turnover_wide %>%
  filter(GeneSampleID %in% rt_samples) %>%
  arrange(match(GeneSampleID, rt_samples))
if (anyDuplicated(turnover_wide$GeneSampleID) > 0) {
  stop("Duplicated GeneSampleID after turnover pivot / GeneSampleID.")
}
turn_mat <- turnover_wide %>%
  select(Turnover_RS, Turnover_BS) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
  as.matrix()
storage.mode(turn_mat) <- "double"
rownames(turn_mat) <- turnover_wide$GeneSampleID
if (any(!is.finite(turn_mat))) {
  stop("Non-finite turnover values detected / .")
}
turn_scaled <- apply(turn_mat, 2, zpop)
turn_scaled <- as.matrix(turn_scaled)
storage.mode(turn_scaled) <- "double"
rownames(turn_scaled) <- turnover_wide$GeneSampleID
turn_pc <- prcomp(
  turn_scaled,
  center = FALSE,
  scale. = FALSE
)
turn_score <- as.numeric(turn_pc$x[, 1])
if (cor(turn_score, rowMeans(turn_scaled), use = "complete.obs") < 0) {
  turn_score <- -turn_score
}
turnover_df <- tibble(
  GeneSampleID = as.character(turnover_wide$GeneSampleID),
  Turnover = turn_score,
  Turnover_PC1_variance = summary(turn_pc)$importance[2, 1]
) %>%
  arrange(match(GeneSampleID, rt_samples))
cat(
  "Turnover PC1 variance = ",
  round(100 * unique(turnover_df$Turnover_PC1_variance), 1),
  "%\n",
  sep = ""
)
# ----------------------------------------------------------------------
# 5. Build the root metabolic module from the 8 pre-selected metabolites
# ----------------------------------------------------------------------
met_ids_high <- c(
  "Com_1057_pos",  # Sakuranetin
  "Com_659_pos",   # Hypoxanthine
  "Com_384_pos",   # Asp-Glu
  "Com_745_pos"    # L-Kynurenine
)
met_ids_low <- c(
  "Com_701_neg",   # Uridine
  "Com_212_neg",   # Catalpol
  "Com_697_neg",   # Tyrosol
  "Com_146_neg"    # 8-O-Acetylharpagide
)
met_ids <- c(met_ids_high, met_ids_low)
met8 <- metab %>%
  filter(FeatureID %in% met_ids)
if (nrow(met8) != 8) {
  stop(
    "Cannot find all eight metabolite features. Found IDs: ",
    paste(met8$FeatureID, collapse = ", ")
  )
}
met_mat <- met8 %>%
  select(FeatureID, all_of(rt_samples)) %>%
  column_to_rownames("FeatureID") %>%
  as.matrix()
# samples x metabolites
met_mat <- t(met_mat)
met_log <- log1p(met_mat)
met_z <- apply(met_log, 2, zpop)
met_z <- as.matrix(met_z)
storage.mode(met_z) <- "double"
# apply() can drop sample dimnames in some situations; restore them explicitly.
# apply(),.
rownames(met_z) <- rt_samples
colnames(met_z) <- colnames(met_log)
# Reverse-code low-yield-associated metabolites
met_z[, met_ids_low] <- -met_z[, met_ids_low]
met_pc <- prcomp(
  met_z,
  center = FALSE,
  scale. = FALSE
)
met_score <- met_pc$x[, 1]
# Orient PC1 so higher scores mean a more high-yield-associated metabolic state.
if (cor(met_score, rowMeans(met_z)) < 0) {
  met_score <- -met_score
  met_pc$rotation[, 1] <- -met_pc$rotation[, 1]
}
metabolic_df <- tibble(
  GeneSampleID = as.character(rt_samples),
  MetabolicModule = as.numeric(met_score),
  MetabolicModule_PC1_variance = summary(met_pc)$importance[2, 1]
) %>%
  arrange(match(GeneSampleID, rt_samples))
met_loading <- tibble(
  FeatureID = rownames(met_pc$rotation),
  PC1_loading = met_pc$rotation[, 1]
)
write_csv(
  met_loading,
  file.path(OUTDIR, "03_metabolic_module_PC1_loadings.csv")
)
cat(
  "Root metabolic module PC1 variance = ",
  round(100 * unique(metabolic_df$MetabolicModule_PC1_variance), 1),
  "%\n",
  sep = ""
)
# ----------------------------------------------------------------------
# 6. Plant traits and final 39-sample analysis table
# ----------------------------------------------------------------------
required_plant_cols <- c(
  "GeneSampleID",
  "Yield_t_ha",
  "RootActivity_ugTTC_gFW_h"
)
missing_plant_cols <- setdiff(required_plant_cols, names(plant))
if (length(missing_plant_cols) > 0) {
  stop(
    "Plant table is missing required columns:\n",
    paste(missing_plant_cols, collapse = ", ")
  )
}
plant_df <- plant %>%
  transmute(
    GeneSampleID = as.character(GeneSampleID),
    Yield = suppressWarnings(as.numeric(`Yield_t_ha`)),
    RootActivity = suppressWarnings(as.numeric(`RootActivity_ugTTC_gFW_h`))
  ) %>%
  filter(GeneSampleID %in% rt_samples) %>%
  arrange(match(GeneSampleID, rt_samples))
if (anyDuplicated(plant_df$GeneSampleID) > 0) {
  stop("Duplicated GeneSampleID in plant table / GeneSampleID.")
}
missing_plant_ids <- setdiff(rt_samples, plant_df$GeneSampleID)
if (length(missing_plant_ids) > 0) {
  stop(
    "Plant table is missing RT samples:\n",
    paste(missing_plant_ids, collapse = ", ")
  )
}
# Explicitly validate every right-hand join table before merging.
# , left_join() .
check_join_table <- function(x, object_name) {
  if (!is.data.frame(x)) {
    stop(object_name, " is not a data.frame/tibble / .")
  }
  if (!"GeneSampleID" %in% names(x)) {
    stop(
      object_name,
      " has no GeneSampleID column / GeneSampleID.
Columns / : ",
      paste(names(x), collapse = ", ")
    )
  }
  x$GeneSampleID <- as.character(x$GeneSampleID)
  if (anyDuplicated(x$GeneSampleID) > 0) {
    dup <- unique(x$GeneSampleID[duplicated(x$GeneSampleID)])
    stop(
      object_name,
      " has duplicated GeneSampleID / GeneSampleID: ",
      paste(dup, collapse = ", ")
    )
  }
  missing_ids <- setdiff(rt_samples, x$GeneSampleID)
  if (length(missing_ids) > 0) {
    stop(
      object_name,
      " is missing metadata samples: ",
      paste(missing_ids, collapse = ", ")
    )
  }
  x %>%
    filter(GeneSampleID %in% rt_samples) %>%
    arrange(match(GeneSampleID, rt_samples))
}
bacteria_df  <- check_join_table(bacteria_df,  "bacteria_df")
turnover_df  <- check_join_table(turnover_df,  "turnover_df")
M1           <- check_join_table(M1,           "M1")
metabolic_df <- check_join_table(metabolic_df, "metabolic_df")
plant_df     <- check_join_table(plant_df,     "plant_df")
cat("\nJoin-table check:\n")
cat("bacteria_df : ", nrow(bacteria_df),  " rows; key = GeneSampleID
", sep = "")
cat("turnover_df : ", nrow(turnover_df),  " rows; key = GeneSampleID
", sep = "")
cat("M1          : ", nrow(M1),           " rows; key = GeneSampleID
", sep = "")
cat("metabolic_df: ", nrow(metabolic_df), " rows; key = GeneSampleID
", sep = "")
cat("plant_df    : ", nrow(plant_df),     " rows; key = GeneSampleID
", sep = "")
dat <- meta %>%
  select(
    GeneSampleID,
    TreatmentCode,
    Treatment,
    Year,
    Replicate
  ) %>%
  mutate(
    GeneSampleID = as.character(GeneSampleID),
    Year = suppressWarnings(as.numeric(Year)),
    Mowing = ifelse(TreatmentCode == "H", 1, 0)
  )
dat <- dplyr::left_join(dat, bacteria_df,  by = "GeneSampleID")
dat <- dplyr::left_join(dat, turnover_df,  by = "GeneSampleID")
dat <- dplyr::left_join(dat, M1,           by = "GeneSampleID")
dat <- dplyr::left_join(dat, metabolic_df, by = "GeneSampleID")
dat <- dplyr::left_join(dat, plant_df,     by = "GeneSampleID")
cat("Final merged rows /  = ", nrow(dat), "
", sep = "")
if (any(!complete.cases(
  dat[, c(
    "BacterialBalance","Turnover","M1",
    "MetabolicModule","RootActivity","Yield"
  )]
))) {
  stop("Missing values after merging the raw files.")
}
dat <- dat %>%
  mutate(
    Year_c = Year - mean(Year),
    Year2_c = Year_c^2,
    BacterialBalance_s = zpop(BacterialBalance),
    Turnover_s = zpop(Turnover),
    M1_s = zpop(M1),
    MetabolicModule_s = zpop(MetabolicModule),
    RootActivity_s = zpop(RootActivity),
    Yield_s = zpop(Yield)
  )
write_csv(
  dat,
  file.path(OUTDIR, "04_SEM_analysis_data_39samples.csv")
)
# ----------------------------------------------------------------------
# 7. Model A and Model B
# ----------------------------------------------------------------------
fit_models <- function(with_bacteria_to_root = FALSE) {
  models <- list(
    Bacteria = lm(
      BacterialBalance_s ~ Mowing + Year_c + Year2_c,
      data = dat
    ),
    Turnover = lm(
      Turnover_s ~ BacterialBalance_s + Year_c + Year2_c,
      data = dat
    ),
    M1 = lm(
      M1_s ~ Turnover_s + BacterialBalance_s + Year_c + Year2_c,
      data = dat
    ),
    Metabolism = lm(
      MetabolicModule_s ~ M1_s + BacterialBalance_s + Year_c + Year2_c,
      data = dat
    ),
    RootActivity = if (!with_bacteria_to_root) {
      lm(
        RootActivity_s ~ MetabolicModule_s + Year_c + Year2_c,
        data = dat
      )
    } else {
      lm(
        RootActivity_s ~ MetabolicModule_s + BacterialBalance_s +
          Year_c + Year2_c,
        data = dat
      )
    },
 # Yield equation
 # Direct effects included in BOTH models
 # Root activity + M1 + RT bacterial state + root metabolism
 # + microbial turnover + mowing -> Yield
 # + M1 + RT +
 # + + ->
    Yield = lm(
      Yield_s ~ RootActivity_s + M1_s +
        BacterialBalance_s + MetabolicModule_s +
        Turnover_s + Mowing +
        Year_c + Year2_c,
      data = dat
    )
  )
  models
}
modA <- fit_models(FALSE)
modB <- fit_models(TRUE)
# ----------------------------------------------------------------------
# 8. Extract path coefficients
# ----------------------------------------------------------------------
extract_paths <- function(models, with_branch = FALSE) {
  rows <- list(
    tibble(
      Source = "Mowing",
      Target = "RT bacterial state",
      Term = "Mowing",
      Model = "Bacteria"
    ),
    tibble(
      Source = "RT bacterial state",
      Target = "Microbial turnover",
      Term = "BacterialBalance_s",
      Model = "Turnover"
    ),
    tibble(
      Source = "RT bacterial state",
      Target = "M1 functional module",
      Term = "BacterialBalance_s",
      Model = "M1"
    ),
    tibble(
      Source = "Microbial turnover",
      Target = "M1 functional module",
      Term = "Turnover_s",
      Model = "M1"
    ),
    tibble(
      Source = "RT bacterial state",
      Target = "Root metabolic module",
      Term = "BacterialBalance_s",
      Model = "Metabolism"
    ),
    tibble(
      Source = "M1 functional module",
      Target = "Root metabolic module",
      Term = "M1_s",
      Model = "Metabolism"
    ),
    tibble(
      Source = "Root metabolic module",
      Target = "Root activity",
      Term = "MetabolicModule_s",
      Model = "RootActivity"
    ),
    tibble(
      Source = "Root activity",
      Target = "Yield",
      Term = "RootActivity_s",
      Model = "Yield"
    ),
    tibble(
      Source = "M1 functional module",
      Target = "Yield",
      Term = "M1_s",
      Model = "Yield"
    ),
    tibble(
      Source = "RT bacterial state",
      Target = "Yield",
      Term = "BacterialBalance_s",
      Model = "Yield"
    ),
    tibble(
      Source = "Root metabolic module",
      Target = "Yield",
      Term = "MetabolicModule_s",
      Model = "Yield"
    ),
    tibble(
      Source = "Microbial turnover",
      Target = "Yield",
      Term = "Turnover_s",
      Model = "Yield"
    ),
    tibble(
      Source = "Mowing",
      Target = "Yield",
      Term = "Mowing",
      Model = "Yield"
    )
  )
  if (with_branch) {
    rows[[length(rows) + 1]] <- tibble(
      Source = "RT bacterial state",
      Target = "Root activity",
      Term = "BacterialBalance_s",
      Model = "RootActivity"
    )
  }
  spec <- bind_rows(rows)
  out <- spec %>%
    rowwise() %>%
    mutate(
      tmp = list(
        get_hc3(
          models[[Model]],
          Term
        )
      ),
      Beta = tmp$Estimate,
      SE_HC3 = tmp$SE_HC3,
      P = tmp$P,
      R2_target = summary(models[[Model]])$r.squared
    ) %>%
    ungroup() %>%
    select(
      Source, Target,
      Beta, SE_HC3, P, R2_target
    )
  out
}
pathsA <- extract_paths(modA, FALSE)
pathsB <- extract_paths(modB, TRUE)
write_csv(
  pathsA,
  file.path(OUTDIR, "05_ModelA_paths.csv")
)
write_csv(
  pathsB,
  file.path(OUTDIR, "06_ModelB_paths.csv")
)
# ----------------------------------------------------------------------
# 9. Fisher's C global fit
# Report both standard OLS d-separation and HC3 sensitivity versions.
# ----------------------------------------------------------------------
basis_formulas <- function(with_branch = FALSE) {
 # IMPORTANT
 # Yield now receives direct paths from Mowing, RT bacterial state,
 # Microbial turnover, M1, Root metabolic module, and Root activity.
 # ,RT,,M1,
 # .
 # Therefore, Mowing_Yield, Turnover_Yield, Bacteria_Yield and
 # Metabolism_Yield are NOT conditional-independence claims anymore
 # and must not enter Fisher's C.
 # ,Fisher's C.
  if (!with_branch) {
    return(list(
      Mowing_Turnover =
        Turnover_s ~ Mowing + BacterialBalance_s + Year_c + Year2_c,
      Mowing_M1 =
        M1_s ~ Mowing + Turnover_s + BacterialBalance_s + Year_c + Year2_c,
      Mowing_Metabolism =
        MetabolicModule_s ~ Mowing + M1_s + BacterialBalance_s +
          Year_c + Year2_c,
      Mowing_Root =
        RootActivity_s ~ Mowing + MetabolicModule_s + Year_c + Year2_c,
      Turnover_Metabolism =
        MetabolicModule_s ~ Turnover_s + M1_s + BacterialBalance_s +
          Year_c + Year2_c,
      Turnover_Root =
        RootActivity_s ~ Turnover_s + MetabolicModule_s +
          Year_c + Year2_c,
      Bacteria_Root =
        RootActivity_s ~ BacterialBalance_s + MetabolicModule_s +
          Year_c + Year2_c,
      M1_Root =
        RootActivity_s ~ M1_s + MetabolicModule_s +
          Year_c + Year2_c
    ))
  }
  list(
    Mowing_Turnover =
      Turnover_s ~ Mowing + BacterialBalance_s + Year_c + Year2_c,
    Mowing_M1 =
      M1_s ~ Mowing + Turnover_s + BacterialBalance_s + Year_c + Year2_c,
    Mowing_Metabolism =
      MetabolicModule_s ~ Mowing + M1_s + BacterialBalance_s +
        Year_c + Year2_c,
    Mowing_Root =
      RootActivity_s ~ Mowing + MetabolicModule_s + BacterialBalance_s +
        Year_c + Year2_c,
    Turnover_Metabolism =
      MetabolicModule_s ~ Turnover_s + M1_s + BacterialBalance_s +
        Year_c + Year2_c,
    Turnover_Root =
      RootActivity_s ~ Turnover_s + MetabolicModule_s +
        BacterialBalance_s + Year_c + Year2_c,
    M1_Root =
      RootActivity_s ~ M1_s + MetabolicModule_s + BacterialBalance_s +
        Year_c + Year2_c
  )
}
fisher_C <- function(formulas, robust = FALSE) {
  ps <- numeric(length(formulas))
  names(ps) <- names(formulas)
  for (i in seq_along(formulas)) {
    f <- formulas[[i]]
    pred <- all.vars(f[[3]])[1]
    m <- lm(f, data = dat)
    if (!robust) {
      ps[i] <- summary(m)$coefficients[pred, "Pr(>|t|)"]
    } else {
      ct <- hc3_table(m)
      ps[i] <- ct$P[ct$Term == pred]
    }
  }
  C <- -2 * sum(log(ps))
  df_C <- 2 * length(ps)
  P_C <- pchisq(C, df = df_C, lower.tail = FALSE)
  tibble(
    Fisher_C = C,
    df = df_C,
    P = P_C
  )
}
fitA_std <- fisher_C(basis_formulas(FALSE), robust = FALSE)
fitA_hc3 <- fisher_C(basis_formulas(FALSE), robust = TRUE)
fitB_std <- fisher_C(basis_formulas(TRUE), robust = FALSE)
fitB_hc3 <- fisher_C(basis_formulas(TRUE), robust = TRUE)
AIC_A <- sum(vapply(modA, AIC, numeric(1)))
AIC_B <- sum(vapply(modB, AIC, numeric(1)))
fit_table <- bind_rows(
  fitA_std %>%
    mutate(
      Model = "A_no_RT_bacteria_to_root_activity",
      Inference = "Standard OLS d-separation",
      Component_AIC_sum = AIC_A
    ),
  fitA_hc3 %>%
    mutate(
      Model = "A_no_RT_bacteria_to_root_activity",
      Inference = "HC3 sensitivity",
      Component_AIC_sum = AIC_A
    ),
  fitB_std %>%
    mutate(
      Model = "B_with_RT_bacteria_to_root_activity",
      Inference = "Standard OLS d-separation",
      Component_AIC_sum = AIC_B
    ),
  fitB_hc3 %>%
    mutate(
      Model = "B_with_RT_bacteria_to_root_activity",
      Inference = "HC3 sensitivity",
      Component_AIC_sum = AIC_B
    )
)
write_csv(
  fit_table,
  file.path(OUTDIR, "07_global_fit_comparison.csv")
)
# ----------------------------------------------------------------------
# 10. Key comparison: attenuation of metabolism -> root activity
# ----------------------------------------------------------------------
compA <- get_hc3(modA$RootActivity, "MetabolicModule_s")
compB_met <- get_hc3(modB$RootActivity, "MetabolicModule_s")
compB_bac <- get_hc3(modB$RootActivity, "BacterialBalance_s")
comparison <- tibble(
  Path = c(
    "Root metabolic module -> Root activity",
    "RT bacterial state -> Root activity"
  ),
  Model_A_Beta = c(compA$Estimate, NA_real_),
  Model_A_P = c(compA$P, NA_real_),
  Model_B_Beta = c(compB_met$Estimate, compB_bac$Estimate),
  Model_B_P = c(compB_met$P, compB_bac$P)
)
write_csv(
  comparison,
  file.path(OUTDIR, "08_key_path_comparison.csv")
)
# Direct effects on Yield in both models
yield_terms <- tibble(
  Path = c(
    "Mowing -> Yield",
    "RT bacterial state -> Yield",
    "Microbial turnover -> Yield",
    "M1 functional module -> Yield",
    "Root metabolic module -> Yield",
    "Root activity -> Yield"
  ),
  Term = c(
    "Mowing",
    "BacterialBalance_s",
    "Turnover_s",
    "M1_s",
    "MetabolicModule_s",
    "RootActivity_s"
  )
)
extract_yield_direct <- function(model, model_name) {
  tab <- hc3_table(model$Yield)
  yield_terms %>%
    left_join(tab, by = "Term") %>%
    transmute(
      Model = model_name,
      Path,
      Beta = Estimate,
      SE_HC3,
      P,
      Yield_R2 = summary(model$Yield)$r.squared,
      Yield_adj_R2 = summary(model$Yield)$adj.r.squared
    )
}
yield_direct_table <- bind_rows(
  extract_yield_direct(modA, "Model A"),
  extract_yield_direct(modB, "Model B")
)
write_csv(
  yield_direct_table,
  file.path(OUTDIR, "08b_direct_effects_on_yield_both_models.csv")
)
# ----------------------------------------------------------------------
# 11. SEM plotting
# ----------------------------------------------------------------------
node_df <- tribble(
  ~Node, ~x, ~y, ~Label,
  "Mowing", 0.06, 0.56, "Mowing",
  "RT bacterial state", 0.26, 0.56, "RT bacterial state",
  "Microbial turnover", 0.47, 0.79, "Microbial turnover",
  "M1 functional module", 0.47, 0.56, "M1 functional module\n(central hub)",
  "Root metabolic module", 0.47, 0.30, "Root metabolic module",
  "Root activity", 0.72, 0.30, "Root activity",
  "Yield", 0.91, 0.56, "Yield"
)
curve_map <- tribble(
  ~Source, ~Target, ~curvature,
  "Mowing", "RT bacterial state", 0,
  "RT bacterial state", "Microbial turnover", -0.25,
  "RT bacterial state", "M1 functional module", 0,
  "Microbial turnover", "M1 functional module", 0,
  "RT bacterial state", "Root metabolic module", 0.20,
  "M1 functional module", "Root metabolic module", 0,
  "RT bacterial state", "Root activity", 0.34,
  "Root metabolic module", "Root activity", 0,
  "Root activity", "Yield", -0.22,
  "M1 functional module", "Yield", 0,
  "RT bacterial state", "Yield", -0.24,
  "Root metabolic module", "Yield", -0.12,
  "Microbial turnover", "Yield", 0.20,
  "Mowing", "Yield", -0.34
)
p_label <- function(p) {
 # Vectorized significance labels for use inside dplyr::mutate().
 # :P.
  dplyr::case_when(
    is.na(p)   ~ "NA",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    p < 0.10   ~ "†",
    TRUE       ~ "n.s."
  )
}
make_sem_plot <- function(path_df, model_name, fit_row) {
  e <- path_df %>%
    left_join(
      node_df %>% select(Source = Node, x1 = x, y1 = y),
      by = "Source"
    ) %>%
    left_join(
      node_df %>% select(Target = Node, x2 = x, y2 = y),
      by = "Target"
    ) %>%
    left_join(curve_map, by = c("Source","Target")) %>%
    mutate(
      curvature = ifelse(is.na(curvature), 0, curvature),
      line_type = case_when(
        P < 0.05 ~ "solid",
        P < 0.10 ~ "dashed",
        TRUE ~ "dotted"
      ),
      label = paste0(
        "β=", sprintf("%.3f", Beta),
        p_label(P)
      ),
      lx = (x1 + x2)/2,
      ly = (y1 + y2)/2
    )
  p <- ggplot()
 # geom_curve() requires one scalar curvature per layer.
 # geom_curve()curvature,.
  for (cv in unique(e$curvature)) {
    ee <- e[abs(e$curvature - cv) < 1e-12, , drop = FALSE]
    p <- p +
      geom_curve(
        data = ee,
        aes(
          x = x1, y = y1,
          xend = x2, yend = y2,
          linetype = line_type
        ),
        curvature = cv,
        linewidth = 0.8,
        arrow = grid::arrow(length = grid::unit(2.8, "mm"), type = "closed"),
        inherit.aes = FALSE
      )
  }
  p <- p +
    geom_label(
      data = node_df,
      aes(x = x, y = y, label = Label),
      size = 3.6,
      label.size = 0.35,
      label.padding = grid::unit(0.35, "lines"),
      fill = "white"
    ) +
    geom_label(
      data = e,
      aes(x = lx, y = ly, label = label),
      size = 2.8,
      label.size = 0,
      fill = "white"
    ) +
    scale_linetype_identity() +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0.13, 0.92),
      clip = "off"
    ) +
    labs(
      title = model_name,
      subtitle = paste0(
        "Fisher's C = ",
        sprintf("%.2f", fit_row$Fisher_C),
        ", df = ", fit_row$df,
        ", P = ", sprintf("%.3f", fit_row$P),
        "; component AIC sum = ",
        sprintf("%.2f", fit_row$Component_AIC_sum)
      )
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9),
      plot.margin = margin(10, 30, 10, 30)
    )
  p
}
fitA_forplot <- fit_table %>%
  filter(
    Model == "A_no_RT_bacteria_to_root_activity",
    Inference == "HC3 sensitivity"
  )
fitB_forplot <- fit_table %>%
  filter(
    Model == "B_with_RT_bacteria_to_root_activity",
    Inference == "HC3 sensitivity"
  )
pA <- make_sem_plot(
  pathsA,
  "Model A: without RT bacterial state -> Root activity",
  fitA_forplot
)
pB <- make_sem_plot(
  pathsB,
  "Model B: with RT bacterial state -> Root activity",
  fitB_forplot
)
ggsave(
  file.path(OUTDIR, "FigS3A_extended_ModelA.png"),
  pA,
  width = 11,
  height = 7,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(OUTDIR, "FigS3A_extended_ModelA.svg"),
  pA,
  width = 11,
  height = 7,
  bg = "white"
)
ggsave(
  file.path(OUTDIR, "FigS3B_extended_ModelB.png"),
  pB,
  width = 11,
  height = 7,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(OUTDIR, "FigS3B_extended_ModelB.svg"),
  pB,
  width = 11,
  height = 7,
  bg = "white"
)
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Final reproducibility record
# ----------------------------------------------------------------------
capture.output(
  sessionInfo(),
  file = file.path(OUTDIR, "sessionInfo_FigS3.txt")
)
cat("\nFigure S3 extended-yield-path sensitivity analysis completed.\n")
