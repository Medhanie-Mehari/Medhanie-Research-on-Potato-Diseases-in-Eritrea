#!/usr/bin/env Rscript
#==============================================================================
# 08_virome_analysis.R
# Target-virus counts, relative abundance, alpha diversity, Bray-Curtis
# dissimilarity and UPGMA clustering from a Kraken2 species count matrix
#
# Authors : <AUTHOR NAMES>
# Paper   : <CITATION / DOI OF PUBLISHED ARTICLE>
# Licence : <LICENCE, e.g. MIT>
#
# Called by 08_virome_pipeline.sh (which sets OUTDIR and METADATA), or run
# directly after setting OUTDIR below.
#
# Input : <OUTDIR>/02_kraken/species_abundance_counts_raw.tsv
# Output: <OUTDIR>/03_target_viruses/  CSV tables T01-T06, figures F01-F03 (PNG + PDF)
#         <OUTDIR>/03_catalogue/       optional detection-only catalogue
#==============================================================================

# -------------------- 0. SETTINGS --------------------
OUTDIR   <- Sys.getenv("OUTDIR",   unset = "<PATH_TO_OUTPUT_FOLDER>")
METADATA <- Sys.getenv("METADATA", unset = "")     # optional sample metadata TSV

# Target viruses: label = regular expression matched (case-insensitive) against
# Kraken2 species names. Include both the current binomial species name and the
# common virus name so that older and newer taxonomy releases both match.
TARGET_VIRUSES <- c(
  "Potato virus Y (Potyvirus)"          = "yituberosi|potato virus y",
  "Potato leafroll virus (Polerovirus)" = "polerovirus plrv|potato leafroll",
  "Potato virus X (Potexvirus)"         = "ecspotati|ecsnidii|potato virus x",
  "Alfalfa mosaic virus (Alfamovirus)"  = "alfamovirus amv|alfalfa mosaic"
)

# Optional: other taxa to list as detection-only (counts table, not analysed).
# Regular expression, or "" to skip, e.g. "Orthotospovirus".
CATALOGUE_PATTERN <- ""

# -------------------- 1. PACKAGES --------------------
required_pkgs <- c("vegan", "ggplot2", "dplyr", "tidyr", "readr", "tibble", "scales", "pheatmap", "ape")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0) stop("Missing R packages: ", paste(missing, collapse = ", "),
                              ". Run 08_virome_setup.sh or install them.")
suppressPackageStartupMessages({
  library(vegan); library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(tibble); library(scales); library(pheatmap); library(ape)
})

# -------------------- 2. PATHS --------------------
if (grepl("<", OUTDIR, fixed = TRUE)) stop("Set OUTDIR (or run via 08_virome_pipeline.sh).")
DIR_CORE <- file.path(OUTDIR, "03_target_viruses")
DIR_CAT  <- file.path(OUTDIR, "03_catalogue")
dir.create(DIR_CORE, showWarnings = FALSE, recursive = TRUE)

raw_file <- file.path(OUTDIR, "02_kraken", "species_abundance_counts_raw.tsv")
if (!file.exists(raw_file)) stop("Matrix file not found: ", raw_file)

# -------------------- 3. LOAD MATRIX, SAMPLE LABELS AND ORDER --------------------
message("Loading species abundance matrix...")
mat_raw  <- readr::read_tsv(raw_file, show_col_types = FALSE)
taxa_all <- as.character(mat_raw[[1]])
counts   <- as.matrix(mat_raw[, -1, drop = FALSE])
rownames(counts) <- taxa_all
storage.mode(counts) <- "numeric"
sample_ids <- colnames(counts)

# Labels, groups and order from metadata (optional)
labels <- setNames(sample_ids, sample_ids)
groups <- setNames(rep(NA_character_, length(sample_ids)), sample_ids)
order_ids <- sort(sample_ids)
if (nzchar(METADATA)) {
  meta <- read.delim(METADATA, check.names = FALSE, colClasses = "character")
  meta <- meta[!startsWith(meta[[1]], "#"), , drop = FALSE]
  ids  <- meta[[1]]
  if ("Label" %in% names(meta)) labels[ids[ids %in% sample_ids]] <- meta$Label[ids %in% sample_ids]
  if ("Group" %in% names(meta)) groups[ids[ids %in% sample_ids]] <- meta$Group[ids %in% sample_ids]
  order_ids <- c(ids[ids %in% sample_ids], setdiff(sort(sample_ids), ids))
  not_found <- setdiff(ids, sample_ids)
  if (length(not_found)) message("Note: metadata samples without a Kraken2 report: ",
                                 paste(not_found, collapse = ", "))
}
counts <- counts[, order_ids, drop = FALSE]
colnames(counts) <- unname(labels[order_ids])
groups <- setNames(unname(groups[order_ids]), colnames(counts))
if (anyDuplicated(colnames(counts))) stop("Sample labels must be unique.")

# -------------------- 4. OPTIONAL DETECTION-ONLY CATALOGUE --------------------
if (nzchar(CATALOGUE_PATTERN)) {
  cat_idx <- grepl(CATALOGUE_PATTERN, taxa_all, ignore.case = TRUE)
  if (any(cat_idx)) {
    dir.create(DIR_CAT, showWarnings = FALSE, recursive = TRUE)
    cat_df <- tibble::rownames_to_column(as.data.frame(counts[cat_idx, , drop = FALSE]), var = "Taxon")
    readr::write_csv(cat_df, file.path(DIR_CAT, "detection_only_counts.csv"))
    message("Detection-only catalogue written to: ", DIR_CAT)
  }
}

# -------------------- 5. TARGET-VIRUS MATRIX --------------------
agg_mat <- matrix(0, nrow = length(TARGET_VIRUSES), ncol = ncol(counts),
                  dimnames = list(names(TARGET_VIRUSES), colnames(counts)))
for (v in names(TARGET_VIRUSES)) {
  idx <- grepl(TARGET_VIRUSES[[v]], taxa_all, ignore.case = TRUE)
  message(sprintf("  %-40s species rows matched: %d", v, sum(idx)))
  agg_mat[v, ] <- colSums(counts[idx, , drop = FALSE])
}

# -------------------- 6. RELATIVE ABUNDANCE AND ALPHA DIVERSITY --------------------
comm   <- t(agg_mat)                      # samples x target viruses
totals <- rowSums(comm)
rel_abund <- comm / ifelse(totals > 0, totals, 1)

alpha_one <- function(x) {
  if (sum(x) == 0) return(c(richness = 0, shannon = 0, gini_simpson = 0))
  p <- x[x > 0] / sum(x)
  c(richness = sum(x > 0), shannon = round(-sum(p * log(p)), 3), gini_simpson = round(1 - sum(p^2), 3))
}
alpha_mat <- t(apply(comm, 1, alpha_one))

# -------------------- 7. TABLES --------------------
message("Writing tables...")
readr::write_csv(tibble::rownames_to_column(as.data.frame(agg_mat), var = "Virus"),
                 file.path(DIR_CORE, "T01_target_virus_counts.csv"))
readr::write_csv(tibble::rownames_to_column(round(as.data.frame(t(rel_abund * 100)), 2), var = "Virus"),
                 file.path(DIR_CORE, "T02_target_virus_relative_abundance_pct.csv"))

alpha_df <- data.frame(Sample = rownames(comm), Group = unname(groups[rownames(comm)]),
                       Observed_Richness = alpha_mat[, "richness"], Shannon = alpha_mat[, "shannon"],
                       Gini_Simpson = alpha_mat[, "gini_simpson"], Total_Target_Reads = as.numeric(totals),
                       stringsAsFactors = FALSE)
readr::write_csv(alpha_df, file.path(DIR_CORE, "T03_alpha_diversity.csv"))

long_df <- as.data.frame(as.table(comm), stringsAsFactors = FALSE)
colnames(long_df) <- c("Sample", "Virus", "Reads")
long_df <- long_df %>% group_by(Sample) %>%
  mutate(Percent_in_Sample = ifelse(sum(Reads) > 0, round(100 * Reads / sum(Reads), 2), 0)) %>%
  ungroup() %>% arrange(Sample, desc(Reads))
readr::write_csv(long_df, file.path(DIR_CORE, "T04_target_virus_long_format.csv"))

pa_matrix <- as.data.frame(ifelse(comm > 0, 1, 0))
readr::write_csv(tibble::rownames_to_column(pa_matrix, var = "Sample"),
                 file.path(DIR_CORE, "T05_presence_absence.csv"))

# -------------------- 8. HEATMAP OF RELATIVE ABUNDANCE --------------------
message("Figure F01: heatmap...")
color_pal <- colorRampPalette(c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404"))(100)
ph <- pheatmap(t(rel_abund * 100), cluster_rows = FALSE, cluster_cols = FALSE, color = color_pal,
               fontsize = 12, fontsize_row = 12, fontsize_col = 11, angle_col = 45,
               display_numbers = TRUE, number_format = "%.1f", fontsize_number = 10,
               number_color = "black", legend = TRUE,
               main = "Relative abundance of target viruses across samples (%)", silent = TRUE)
fig_w <- max(8, 0.75 * ncol(agg_mat) + 4)
ggsave(file.path(DIR_CORE, "F01_target_virus_heatmap.png"), ph$gtable, width = fig_w, height = 5.5, dpi = 300)
ggsave(file.path(DIR_CORE, "F01_target_virus_heatmap.pdf"), ph$gtable, width = fig_w, height = 5.5)

# -------------------- 9. ALPHA DIVERSITY LINE GRAPH --------------------
message("Figure F02: alpha diversity...")
alpha_plot_df <- alpha_df %>% mutate(Sample = factor(Sample, levels = rownames(comm)))
max_shannon <- max(alpha_plot_df$Shannon, na.rm = TRUE); if (!is.finite(max_shannon) || max_shannon == 0) max_shannon <- 1
max_simpson <- max(alpha_plot_df$Gini_Simpson, na.rm = TRUE); if (!is.finite(max_simpson) || max_simpson == 0) max_simpson <- 1
scale_factor <- max_shannon / max_simpson

p_line <- ggplot(alpha_plot_df, aes(x = Sample, group = 1)) +
  geom_line(aes(y = Shannon, color = "Shannon"), linewidth = 1.3) +
  geom_point(aes(y = Shannon, color = "Shannon"), size = 3.8) +
  geom_line(aes(y = Gini_Simpson * scale_factor, color = "Gini-Simpson"), linewidth = 1.3) +
  geom_point(aes(y = Gini_Simpson * scale_factor, color = "Gini-Simpson"), size = 3.8) +
  scale_y_continuous(name = "Shannon diversity index (H')", limits = c(0, max_shannon * 1.18),
                     sec.axis = sec_axis(~ . / scale_factor, name = "Gini–Simpson index (1 - D)")) +
  scale_color_manual(name = "Diversity index", values = c("Shannon" = "#1F78B4", "Gini-Simpson" = "#E31A1C")) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "black", size = 11),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(color = "#1F78B4", face = "bold", size = 12),
        axis.title.y.right = element_text(color = "#E31A1C", face = "bold", size = 12),
        legend.position = "top", legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 11),
        panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")) +
  labs(title = "Alpha diversity of target-virus communities across samples", x = "Sample")

# Group boundaries (samples are shown in metadata order)
g <- groups[levels(alpha_plot_df$Sample)]
if (any(!is.na(g)) && length(unique(g)) > 1) {
  runs   <- rle(ifelse(is.na(g), "", g))
  ends   <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  for (b in head(ends, -1)) p_line <- p_line +
    geom_vline(xintercept = b + 0.5, linetype = "dotted", color = "grey40", linewidth = 0.8)
  p_line <- p_line + annotate("text", x = (starts + ends) / 2, y = max_shannon * 1.12,
                              label = runs$values, size = 3.8, fontface = "italic", color = "grey30")
}
fig_w <- max(8, 0.7 * nrow(comm) + 2)
ggsave(file.path(DIR_CORE, "F02_alpha_diversity.png"), p_line, width = fig_w, height = 6.2, dpi = 300)
ggsave(file.path(DIR_CORE, "F02_alpha_diversity.pdf"), p_line, width = fig_w, height = 6.2)

# -------------------- 10. BRAY-CURTIS AND UPGMA CLUSTERING --------------------
message("Figure F03: UPGMA clustering...")
has_reads <- totals > 0
if (sum(has_reads) < 3) {
  message("Fewer than 3 samples with target-virus reads: clustering skipped.")
} else {
  if (any(!has_reads)) message("Samples without target-virus reads excluded from clustering: ",
                               paste(rownames(comm)[!has_reads], collapse = ", "))
  bc <- vegdist(rel_abund[has_reads, , drop = FALSE], method = "bray")
  readr::write_csv(tibble::rownames_to_column(as.data.frame(as.matrix(bc)), var = "Sample"),
                   file.path(DIR_CORE, "T06_bray_curtis_matrix.csv"))
  hc  <- hclust(bc, method = "average")
  phy <- as.phylo(hc)
  phy$edge.length <- 2 * phy$edge.length   # as.phylo halves hclust heights; restore dissimilarity scale

  grp_levels <- unique(na.omit(groups))
  grp_cols   <- setNames(c("black", "#990000", "#1F78B4", "#33A02C", "#6A3D9A")[seq_along(grp_levels)], grp_levels)
  tip_cols   <- ifelse(is.na(groups[phy$tip.label]), "black", grp_cols[groups[phy$tip.label]])

  draw_tree <- function() {
    par(mar = c(5.5, 3, 4, 9), family = "sans")
    plot(phy, type = "phylogram", direction = "leftwards", tip.color = tip_cols, cex = 1.05,
         font = 2, edge.width = 1.8,
         main = "UPGMA clustering of samples based on\nBray-Curtis dissimilarity of target viruses")
    axisPhylo(side = 1, backward = TRUE)   # 0 at the tips = merge heights
    title(xlab = "Bray-Curtis dissimilarity", line = 3, cex.lab = 1.1, font.lab = 2)
    if (length(grp_levels) > 1)
      legend("topright", legend = grp_levels, text.col = grp_cols, bty = "n", text.font = 2)
  }
  png(file.path(DIR_CORE, "F03_UPGMA_clustering.png"), width = 3000, height = 2100, res = 300)
  draw_tree(); invisible(dev.off())
  pdf(file.path(DIR_CORE, "F03_UPGMA_clustering.pdf"), width = 10, height = 7)
  draw_tree(); invisible(dev.off())
}

message("=== All tables and figures completed ===")
message("Destination: ", DIR_CORE)
