#!/usr/bin/env Rscript

# ============================================================
# Script: 08_program_scores_cell_states.R
#
# Purpose:
#   Scores selected cell-state marker programs in bulk RNA-seq
#   samples using a manual UP-only singscore-like rank score.
#
#   Marker programs are derived from the spinal cord injury
#   single-cell atlas marker table and collapsed into broader
#   biologically interpretable programs.
#
# Inputs:
#   - results/limma_voom/objects/logcpm_annotated.tsv
#   - results/limma_voom/objects/sample_metadata.tsv
#   - data/references/marker_panels_from_atlas_MOESM4.tsv
#
# Outputs:
#   - results/program_scores/tables/atlas_collapsed_marker_panel.tsv
#   - results/program_scores/tables/singscore_collapsed_program_scores.tsv
#   - results/program_scores/tables/program_scores_limma_results.tsv
#   - results/program_scores/tables/program_group_delta_stats.tsv
#   - results/program_scores/tables/program_interaction_results.tsv
#   - results/program_scores/figures/program_group_mean_heatmap.pdf/.png
#   - results/program_scores/figures/program_delta_dotplots.pdf/.png
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(patchwork)
  library(yaml)
})

# -----------------------------
# Parameters
# -----------------------------

MAX_GENES_PER_PROGRAM <- 40
MIN_GENES_PER_PROGRAM <- 8

GROUP_ORDER <- c(
  "WT_SHAM_d7",
  "WT_SCI_d7",
  "KO_SHAM_d7",
  "KO_SCI_d7",
  "DTR_PBS_d14",
  "DTR_DTX_d14"
)

CONTRASTS <- c(
  "WT_SCI_vs_SHAM_dpi7",
  "KO_SCI_vs_SHAM_dpi7",
  "DTR_DTX_vs_PBS_dpi14",
  "Injury_interaction_WTvsKO"
)

# -----------------------------
# Paths
# -----------------------------

load_paths <- function(config_file = "config/paths.yml") {
  if (!file.exists(config_file)) stop("Missing config/paths.yml.")
  yaml::read_yaml(config_file)
}

as_project_path <- function(path, paths) {
  if (grepl("^/", path)) path else file.path(paths$project$root, path)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

paths <- load_paths()

f_logcpm <- file.path(as_project_path(paths$outputs$limma_voom, paths), "objects", "logcpm_annotated.tsv")
f_meta <- file.path(as_project_path(paths$outputs$limma_voom, paths), "objects", "sample_metadata.tsv")
f_markers <- as_project_path(paths$references$marker_panels, paths)

out_base <- as_project_path(paths$outputs$program_scores, paths)
out_tab <- ensure_dir(file.path(out_base, "tables"))
out_fig <- ensure_dir(file.path(out_base, "figures"))

stopifnot(file.exists(f_logcpm))
stopifnot(file.exists(f_meta))
stopifnot(file.exists(f_markers))

# -----------------------------
# Load expression and metadata
# -----------------------------

logcpm <- data.table::fread(f_logcpm) %>% as.data.frame()
meta <- data.table::fread(f_meta) %>% as.data.frame()

stopifnot(all(c("gene_id", "gene_symbol", "gene_name") %in% colnames(logcpm)))
stopifnot(all(c("sample_id", "group_short") %in% colnames(meta)))
stopifnot(all(meta$sample_id %in% colnames(logcpm)))

sample_cols <- meta$sample_id

expr_df <- logcpm %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_symbol, .keep_all = TRUE)

X <- as.matrix(expr_df[, sample_cols, drop = FALSE])
rownames(X) <- expr_df$gene_symbol
storage.mode(X) <- "double"

# -----------------------------
# Load and collapse atlas marker programs
# -----------------------------

markers_raw <- data.table::fread(f_markers) %>% as.data.frame()
stopifnot(all(c("program", "marker_genes") %in% colnames(markers_raw)))

markers_raw <- markers_raw %>%
  mutate(
    program = as.character(program),
    genes = str_split(marker_genes, "\\s*,\\s*")
  ) %>%
  select(program, genes)

program_map <- list(
  Astrocyte_homeostatic = c("Astrocytes 1", "Astrocytes 2"),
  Astrocyte_reactive = c("Reactive Astrocyte", "WM Astrocytes"),

  Microglia_homeostatic = c("Microglia 1", "Microglia 2"),
  Microglia_activated_A = c("Activated Microglia A"),
  Microglia_activated_B = c("Activated Microglia B"),

  Macrophage = c("Macrophages"),
  NK_T = c("NK/T cells"),

  OPC_COP = c("OPC", "COP"),
  Oligodendrocyte_myelinating = c("NFOL", "MFOL", "MOL-1", "MOL-2"),

  Endothelial = c("Endothelial"),
  Perivascular_border = c("Pericytes", "Leptomeninges"),

  Neuronal_synaptic = c(
    "Cpne4", "Maf", "Reln", "Rreb1", "Sox5", "Megf11",
    "ME", "VE", "Adamts5", "Cdh3", "Pdyn", "Npy", "Chat",
    "MI", "VI", "MN"
  )
)

build_program <- function(source_programs, markers_df, max_genes = MAX_GENES_PER_PROGRAM) {
  sub <- markers_df %>%
    filter(program %in% source_programs)

  if (nrow(sub) == 0) return(character(0))

  genes_all <- unlist(sub$genes)

  gene_df <- data.frame(
    gene_symbol = names(table(genes_all)),
    recurrence = as.integer(table(genes_all)),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(recurrence), gene_symbol)

  genes <- head(gene_df$gene_symbol, max_genes)

  genes <- genes[!grepl("^mt-", genes, ignore.case = TRUE)]
  genes <- genes[!grepl("^Rpl|^Rps", genes)]

  unique(genes)
}

collapsed_markers <- bind_rows(lapply(names(program_map), function(program_name) {
  genes <- build_program(program_map[[program_name]], markers_raw)

  data.frame(
    program = program_name,
    gene_symbol = genes,
    stringsAsFactors = FALSE
  )
}))

collapsed_markers <- collapsed_markers %>%
  filter(gene_symbol %in% rownames(X)) %>%
  group_by(program) %>%
  filter(n() >= MIN_GENES_PER_PROGRAM) %>%
  ungroup()

marker_panel <- collapsed_markers %>%
  group_by(program) %>%
  summarise(
    marker_genes = paste(gene_symbol, collapse = ","),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  arrange(program)

readr::write_tsv(
  marker_panel,
  file.path(out_tab, "atlas_collapsed_marker_panel.tsv")
)

# -----------------------------
# Manual UP-only singscore-like scoring
# -----------------------------

markers_list <- collapsed_markers %>%
  group_by(program) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  arrange(program)

R <- apply(X, 2, rank, ties.method = "average")
rownames(R) <- rownames(X)
N <- nrow(R)

score_up <- function(rank_matrix, genes) {
  genes <- intersect(genes, rownames(rank_matrix))
  k <- length(genes)

  if (k < 5) return(rep(NA_real_, ncol(rank_matrix)))
  if (N == k) return(rep(NA_real_, ncol(rank_matrix)))

  mu <- colMeans(rank_matrix[genes, , drop = FALSE], na.rm = TRUE)
  2 * (mu - (N + 1) / 2) / (N - k)
}

score_mat <- t(sapply(seq_len(nrow(markers_list)), function(i) {
  score_up(R, markers_list$genes[[i]])
}))

rownames(score_mat) <- markers_list$program
colnames(score_mat) <- colnames(X)
score_mat <- score_mat[, meta$sample_id, drop = FALSE]

stopifnot(!any(is.na(score_mat)))

readr::write_tsv(
  data.frame(program = rownames(score_mat), score_mat, check.names = FALSE),
  file.path(out_tab, "singscore_collapsed_program_scores.tsv")
)

# -----------------------------
# limma testing on program scores
# -----------------------------

meta$group_short <- factor(meta$group_short, levels = GROUP_ORDER)

design <- model.matrix(~ 0 + group_short, data = meta)
colnames(design) <- levels(meta$group_short)

contrast_matrix <- limma::makeContrasts(
  WT_SCI_vs_SHAM_dpi7 = WT_SCI_d7 - WT_SHAM_d7,
  KO_SCI_vs_SHAM_dpi7 = KO_SCI_d7 - KO_SHAM_d7,
  DTR_DTX_vs_PBS_dpi14 = DTR_DTX_d14 - DTR_PBS_d14,
  Injury_interaction_WTvsKO =
    (KO_SCI_d7 - KO_SHAM_d7) - (WT_SCI_d7 - WT_SHAM_d7),
  levels = design
)

stopifnot(all(CONTRASTS %in% colnames(contrast_matrix)))

fit <- limma::lmFit(score_mat, design)
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2)

program_results <- bind_rows(lapply(CONTRASTS, function(contrast_name) {
  tt <- limma::topTable(
    fit2,
    coef = contrast_name,
    number = Inf,
    sort.by = "none"
  )

  tt$program <- rownames(tt)
  tt$contrast <- contrast_name
  tt
})) %>%
  relocate(contrast, program)

readr::write_tsv(
  program_results,
  file.path(out_tab, "program_scores_limma_results.tsv")
)

interaction_results <- program_results %>%
  filter(contrast == "Injury_interaction_WTvsKO") %>%
  arrange(adj.P.Val) %>%
  mutate(
    direction = ifelse(
      logFC > 0,
      "KO > WT injury effect",
      "WT > KO injury effect"
    ),
    significant_fdr_0_05 = adj.P.Val < 0.05
  )

readr::write_tsv(
  interaction_results,
  file.path(out_tab, "program_interaction_results.tsv")
)

# -----------------------------
# Simple group delta summaries
# -----------------------------

get_group_scores <- function(group_name) {
  samples <- meta$sample_id[meta$group_short == group_name]
  score_mat[, samples, drop = FALSE]
}

contrast_map <- tibble::tribble(
  ~contrast, ~baseline_group, ~comparison_group,
  "WT_SCI_vs_SHAM_dpi7", "WT_SHAM_d7", "WT_SCI_d7",
  "KO_SCI_vs_SHAM_dpi7", "KO_SHAM_d7", "KO_SCI_d7",
  "DTR_DTX_vs_PBS_dpi14", "DTR_PBS_d14", "DTR_DTX_d14"
)

group_delta_stats <- bind_rows(lapply(seq_len(nrow(contrast_map)), function(i) {
  contrast_name <- contrast_map$contrast[i]
  baseline_group <- contrast_map$baseline_group[i]
  comparison_group <- contrast_map$comparison_group[i]

  baseline_mat <- get_group_scores(baseline_group)
  comparison_mat <- get_group_scores(comparison_group)

  data.frame(
    contrast = contrast_name,
    program = rownames(score_mat),
    mean_baseline = rowMeans(baseline_mat, na.rm = TRUE),
    mean_comparison = rowMeans(comparison_mat, na.rm = TRUE),
    delta_comparison_minus_baseline =
      rowMeans(comparison_mat, na.rm = TRUE) - rowMeans(baseline_mat, na.rm = TRUE),
    p_value = sapply(rownames(score_mat), function(program) {
      t.test(comparison_mat[program, ], baseline_mat[program, ])$p.value
    }),
    stringsAsFactors = FALSE
  )
})) %>%
  group_by(contrast) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

readr::write_tsv(
  group_delta_stats,
  file.path(out_tab, "program_group_delta_stats.tsv")
)

# -----------------------------
# Figures
# -----------------------------

group_means <- do.call(cbind, lapply(GROUP_ORDER, function(group_name) {
  samples <- meta$sample_id[meta$group_short == group_name]
  rowMeans(score_mat[, samples, drop = FALSE], na.rm = TRUE)
}))
colnames(group_means) <- GROUP_ORDER

group_means_z <- t(scale(t(group_means)))
group_means_z <- group_means_z[order(apply(group_means_z, 1, sd, na.rm = TRUE), decreasing = TRUE), ]

pheatmap::pheatmap(
  group_means_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA,
  main = "Program activity across groups",
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45,
  filename = file.path(out_fig, "program_group_mean_heatmap.pdf"),
  width = 6.5,
  height = 4.8
)

png(
  filename = file.path(out_fig, "program_group_mean_heatmap.png"),
  width = 6.5,
  height = 4.8,
  units = "in",
  res = 300
)
pheatmap::pheatmap(
  group_means_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA,
  main = "Program activity across groups",
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45
)
dev.off()

plot_delta <- group_delta_stats %>%
  mutate(
    direction = ifelse(delta_comparison_minus_baseline > 0, "Up", "Down"),
    minus_log10_fdr = -log10(p_adj)
  )

make_delta_plot <- function(contrast_name) {
  df <- plot_delta %>%
    filter(contrast == contrast_name) %>%
    arrange(delta_comparison_minus_baseline) %>%
    mutate(program = factor(program, levels = program))

  ggplot(df, aes(x = delta_comparison_minus_baseline, y = program)) +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_point(
      aes(size = minus_log10_fdr, fill = direction),
      shape = 21,
      color = "black",
      stroke = 0.25
    ) +
    scale_fill_manual(values = c(Up = "#C95D63", Down = "#4C78A8")) +
    scale_size(range = c(2, 5)) +
    labs(
      title = contrast_name,
      x = "Δ program score",
      y = NULL
    ) +
    theme_classic(base_size = 9) +
    theme(legend.position = "right")
}

p_delta <- make_delta_plot("WT_SCI_vs_SHAM_dpi7") /
  make_delta_plot("KO_SCI_vs_SHAM_dpi7") /
  make_delta_plot("DTR_DTX_vs_PBS_dpi14")

ggsave(
  file.path(out_fig, "program_delta_dotplots.pdf"),
  p_delta,
  width = 7,
  height = 10
)

ggsave(
  file.path(out_fig, "program_delta_dotplots.png"),
  p_delta,
  width = 7,
  height = 10,
  dpi = 300
)

message("Program score analysis complete.")
message("Tables: ", out_tab)
message("Figures: ", out_fig)
