# ============================================================
# Script: 01_limma_voom_de.R
#
# Purpose:
#   Runs the core bulk RNA-seq differential expression analysis for
#   the Gpnmb spinal cord injury project using edgeR + limma-voom.
#
# Inputs:
#   - config/paths.yml
#   - raw featureCounts gene-by-sample count matrix
#   - raw sample metadata table with columns: name, group1
#
# Outputs:
#   - filtered TMM-normalized logCPM matrix with gene annotations
#   - cleaned sample metadata used by downstream scripts
#   - limma-voom model objects
#   - annotated differential expression tables for all contrasts
#   - minimal run summaries for reproducibility
#
# Notes:
#   This script intentionally keeps plotting minimal. Figure-facing QC,
#   PCA, volcano plots, and concordance plots are produced downstream.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(janitor)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

set.seed(1)

source("R/io_helpers.R")
paths <- load_paths("config/paths.yml")

as_project_path <- function(path) {
  if (grepl("^/", path)) path else file.path(paths$project$root, path)
}

# -----------------------------
# 1. Input and output paths
# -----------------------------

counts_path <- as_project_path(paths$inputs$raw_counts)
metadata_path <- as_project_path(paths$inputs$raw_metadata)

out_dir <- as_project_path(paths$outputs$limma_voom)
dir_qc  <- file.path(out_dir, "qc")
dir_deg <- file.path(out_dir, "deg")
dir_obj <- file.path(out_dir, "objects")

ensure_dir(dir_qc)
ensure_dir(dir_deg)
ensure_dir(dir_obj)

check_file_exists(counts_path, "raw count matrix")
check_file_exists(metadata_path, "raw sample metadata")

# -----------------------------
# 2. Read counts and metadata
# -----------------------------

counts_df <- readr::read_tsv(counts_path, show_col_types = FALSE) |>
  janitor::clean_names()

meta_raw <- readr::read_tsv(metadata_path, show_col_types = FALSE) |>
  janitor::clean_names()

stopifnot("geneid" %in% names(counts_df))
stopifnot(all(c("name", "group1") %in% names(meta_raw)))

count_mat <- counts_df |>
  dplyr::select(geneid, where(is.numeric)) |>
  dplyr::distinct(geneid, .keep_all = TRUE) |>
  tibble::column_to_rownames("geneid") |>
  as.matrix()

storage.mode(count_mat) <- "integer"

meta <- meta_raw |>
  dplyr::transmute(
    sample_raw = as.character(name),
    sample = janitor::make_clean_names(sample_raw),
    group_raw = as.character(group1),
    group = factor(tolower(group_raw))
  )

stopifnot(!anyDuplicated(rownames(count_mat)))
stopifnot(!anyDuplicated(meta$sample))

missing_in_counts <- setdiff(meta$sample, colnames(count_mat))
missing_in_meta <- setdiff(colnames(count_mat), meta$sample)

if (length(missing_in_counts) > 0 || length(missing_in_meta) > 0) {
  stop(
    "Sample name mismatch between metadata and counts.\n",
    "Missing in counts: ", paste(missing_in_counts, collapse = ", "), "\n",
    "Missing in metadata: ", paste(missing_in_meta, collapse = ", ")
  )
}

count_mat <- count_mat[, meta$sample, drop = FALSE]
stopifnot(identical(colnames(count_mat), meta$sample))

# -----------------------------
# 3. Design matrix and contrasts
# -----------------------------

# These group names are the cleaned group labels produced from group1.
# The script stops if the metadata do not contain exactly the expected labels.
expected_groups <- c(
  "wtshamadultdpi7",
  "wtsciadultdpi7",
  "gpnmbkoshamadultdpi7",
  "gpnmbkosciadultdpi7",
  "gdtpbsd7to13adultscidpi14",
  "gdtdtd7to13adultscidpi14"
)

missing_groups <- setdiff(expected_groups, levels(meta$group))
if (length(missing_groups) > 0) {
  stop("Missing expected group labels: ", paste(missing_groups, collapse = ", "))
}

meta$group <- factor(meta$group, levels = expected_groups)

design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

contrasts <- limma::makeContrasts(
  WT_SCI_vs_SHAM_dpi7 =
    wtsciadultdpi7 - wtshamadultdpi7,

  KO_SCI_vs_SHAM_dpi7 =
    gpnmbkosciadultdpi7 - gpnmbkoshamadultdpi7,

  DTR_DTX_vs_PBS_dpi14 =
    gdtdtd7to13adultscidpi14 - gdtpbsd7to13adultscidpi14,

  Injury_interaction_WTvsKO =
    (wtsciadultdpi7 - wtshamadultdpi7) -
    (gpnmbkosciadultdpi7 - gpnmbkoshamadultdpi7),

  # Exploratory cross-cohort/timepoint contrasts. Do not use these
  # as primary causal evidence because cohort and timepoint differ.
  DTR_DTX_vs_WT_SHAM_exploratory =
    gdtdtd7to13adultscidpi14 - wtshamadultdpi7,

  DTR_PBS_vs_WT_SHAM_exploratory =
    gdtpbsd7to13adultscidpi14 - wtshamadultdpi7,

  levels = design
)

readr::write_tsv(as.data.frame(table(meta$group)), file.path(dir_qc, "group_sizes.tsv"))
write.csv(design, file.path(dir_obj, "design_matrix.csv"), row.names = TRUE)
write.csv(contrasts, file.path(dir_obj, "contrast_matrix.csv"), row.names = TRUE)

# -----------------------------
# 4. edgeR filtering, TMM normalization, and limma-voom model
# -----------------------------

y <- edgeR::DGEList(counts = count_mat, group = meta$group)
keep <- edgeR::filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y, method = "TMM")

v <- limma::voomWithQualityWeights(y, design, plot = FALSE)
fit <- limma::lmFit(v, design)
fit2 <- limma::contrasts.fit(fit, contrasts)
fit2 <- limma::eBayes(fit2, robust = TRUE)

filtering_summary <- tibble::tibble(
  n_genes_raw = nrow(count_mat),
  n_genes_kept = nrow(y),
  n_samples = ncol(count_mat)
)
readr::write_tsv(filtering_summary, file.path(dir_obj, "filtering_summary.tsv"))

# -----------------------------
# 5. Gene annotation
# -----------------------------

# Ensembl IDs remain the primary gene identifiers. Symbols/names are added
# only for readability in tables and figures.
gene_id <- rownames(y)
gene_id_stripped <- sub("\\..*$", "", gene_id)

map_df <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(gene_id_stripped),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
) |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename(
    gene_id_stripped = ENSEMBL,
    gene_symbol = SYMBOL,
    gene_name = GENENAME
  )

gene_map <- tibble::tibble(
  gene_id = gene_id,
  gene_id_stripped = gene_id_stripped
) |>
  dplyr::left_join(map_df, by = "gene_id_stripped")

mapping_qc <- tibble::tibble(
  n_genes = nrow(gene_map),
  n_mapped_symbol = sum(!is.na(gene_map$gene_symbol)),
  frac_mapped_symbol = mean(!is.na(gene_map$gene_symbol)),
  org_mm_eg_db_version = as.character(packageVersion("org.Mm.eg.db"))
)

readr::write_tsv(gene_map, file.path(dir_obj, "gene_map_ensembl_to_symbol.tsv"))
readr::write_tsv(mapping_qc, file.path(dir_obj, "gene_mapping_qc.tsv"))

# -----------------------------
# 6. Export normalized expression and model objects
# -----------------------------

logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)

logcpm_annotated <- as.data.frame(logcpm) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::left_join(gene_map, by = "gene_id") |>
  dplyr::select(gene_id, gene_symbol, gene_name, dplyr::everything())

group_labels <- c(
  wtshamadultdpi7 = "WT_SHAM_d7",
  wtsciadultdpi7 = "WT_SCI_d7",
  gpnmbkoshamadultdpi7 = "KO_SHAM_d7",
  gpnmbkosciadultdpi7 = "KO_SCI_d7",
  gdtpbsd7to13adultscidpi14 = "DTR_PBS_d14",
  gdtdtd7to13adultscidpi14 = "DTR_DTX_d14"
)

group_colors <- c(
  WT_SHAM_d7 = "#A6CEE3",
  WT_SCI_d7 = "#1F78B4",
  KO_SHAM_d7 = "#FDBF6F",
  KO_SCI_d7 = "#FF7F00",
  DTR_PBS_d14 = "#B2DF8A",
  DTR_DTX_d14 = "#33A02C"
)

meta_export <- meta |>
  dplyr::mutate(
    group_short = unname(group_labels[as.character(group)]),
    color = unname(group_colors[group_short])
  ) |>
  dplyr::group_by(group_short) |>
  dplyr::mutate(
    sample_index = dplyr::row_number(),
    sample_label = paste0(sample_index, "_", group_short)
  ) |>
  dplyr::ungroup()

stopifnot(!any(is.na(meta_export$group_short)))
stopifnot(!any(is.na(meta_export$color)))

readr::write_tsv(logcpm_annotated, file.path(dir_obj, "logcpm_annotated.tsv"))
readr::write_tsv(meta_export, file.path(dir_obj, "sample_metadata.tsv"))

saveRDS(v, file.path(dir_obj, "voom_object.rds"))
saveRDS(fit2, file.path(dir_obj, "fit2_ebayes.rds"))

# -----------------------------
# 7. Export differential expression tables
# -----------------------------

for (contrast_name in colnames(contrasts)) {
  tt <- limma::topTable(fit2, coef = contrast_name, number = Inf, sort.by = "P") |>
    tibble::rownames_to_column("gene_id") |>
    dplyr::left_join(gene_map, by = "gene_id") |>
    dplyr::select(gene_id, gene_symbol, gene_name, dplyr::everything())

  readr::write_tsv(
    tt,
    file.path(dir_deg, paste0("deg_", contrast_name, "_annotated.tsv"))
  )
}

# -----------------------------
# 8. Reproducibility record
# -----------------------------

run_summary <- list(
  script = "01_limma_voom_de.R",
  counts_path = counts_path,
  metadata_path = metadata_path,
  n_genes_raw = nrow(count_mat),
  n_genes_kept = nrow(y),
  n_samples = ncol(count_mat),
  contrasts = colnames(contrasts),
  voom_quality_weights = TRUE,
  ebayes_robust = TRUE
)

saveRDS(run_summary, file.path(dir_obj, "run_summary.rds"))
capture.output(sessionInfo(), file = file.path(dir_obj, "sessionInfo.txt"))

message("limma-voom analysis complete.")
message("Outputs written to: ", out_dir)
