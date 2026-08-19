suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

counts_path <- snakemake@input[["counts"]]
samples_path <- snakemake@input[["samples"]]
contrasts_path <- snakemake@input[["contrasts"]]
out_dir <- snakemake@output[[1]]
log_path <- snakemake@log[[1]]
design_terms <- unlist(snakemake@params[["design"]])
alpha <- as.numeric(snakemake@params[["alpha"]])
min_total_count <- as.integer(snakemake@params[["min_total_count"]])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_path, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

counts_df <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"gene_id" %in% names(counts_df)) stop("Count matrix must contain gene_id")
rownames(counts_df) <- counts_df$gene_id
counts_df$gene_id <- NULL
count_matrix <- as.matrix(counts_df)
storage.mode(count_matrix) <- "integer"

metadata <- read.delim(samples_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample_id" %in% names(metadata)) stop("Sample sheet must contain sample_id")
if (anyDuplicated(metadata$sample_id)) stop("Duplicate sample_id values in metadata")
rownames(metadata) <- metadata$sample_id

missing_counts <- setdiff(metadata$sample_id, colnames(count_matrix))
extra_counts <- setdiff(colnames(count_matrix), metadata$sample_id)
if (length(missing_counts)) stop("Samples missing from counts: ", paste(missing_counts, collapse = ", "))
if (length(extra_counts)) stop("Count columns missing from metadata: ", paste(extra_counts, collapse = ", "))
count_matrix <- count_matrix[rownames(count_matrix) != "", metadata$sample_id, drop = FALSE]

missing_terms <- setdiff(design_terms, names(metadata))
if (length(missing_terms)) stop("Design terms missing from metadata: ", paste(missing_terms, collapse = ", "))
for (term in design_terms) {
  metadata[[term]] <- droplevels(factor(metadata[[term]]))
  if (nlevels(metadata[[term]]) < 2) stop("Design term has fewer than two levels: ", term)
}

keep <- rowSums(count_matrix) >= min_total_count
if (sum(keep) < 2) stop("Too few genes remain after count filtering")
count_matrix <- count_matrix[keep, , drop = FALSE]
design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = design_formula
)
dds <- DESeq(dds, quiet = TRUE)

normalized <- counts(dds, normalized = TRUE)
write.table(
  data.frame(gene_id = rownames(normalized), normalized, check.names = FALSE),
  file = file.path(out_dir, "normalized_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
write.table(
  data.frame(gene_id = rownames(assay(vsd)), assay(vsd), check.names = FALSE),
  file = file.path(out_dir, "vst_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

pca_group <- tail(design_terms, 1)
pca <- plotPCA(vsd, intgroup = pca_group, returnData = TRUE)
percent_var <- round(100 * attr(pca, "percentVar"))
pca$sample_id <- rownames(pca)
write.table(pca, file.path(out_dir, "pca.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot(pca, aes(x = PC1, y = PC2, color = .data[[pca_group]], label = sample_id)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, check_overlap = TRUE, size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw(base_size = 12) +
  labs(color = pca_group, title = "Sample PCA")
ggsave(file.path(out_dir, "pca.png"), p, width = 8, height = 6, dpi = 180)

sample_dist <- as.matrix(dist(t(assay(vsd))))
png(file.path(out_dir, "sample_distance.png"), width = 1800, height = 1600, res = 200)
pheatmap(sample_dist, main = "Sample distance", border_color = NA)
dev.off()

contrasts <- read.delim(contrasts_path, check.names = FALSE, stringsAsFactors = FALSE)
summary_rows <- list()

for (i in seq_len(nrow(contrasts))) {
  contrast_id <- contrasts$contrast_id[i]
  factor_name <- contrasts$factor[i]
  numerator <- contrasts$numerator[i]
  denominator <- contrasts$denominator[i]
  label <- if ("label" %in% names(contrasts)) contrasts$label[i] else contrast_id

  if (!factor_name %in% design_terms) stop("Contrast factor is not in design: ", factor_name)
  result <- results(
    dds,
    contrast = c(factor_name, numerator, denominator),
    alpha = alpha
  )
  result_df <- data.frame(gene_id = rownames(result), as.data.frame(result), check.names = FALSE)
  result_df <- result_df[order(result_df$padj, na.last = TRUE), ]
  contrast_dir <- file.path(out_dir, contrast_id)
  dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    result_df,
    file.path(contrast_dir, "results.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  significant <- !is.na(result_df$padj) & result_df$padj < alpha
  up <- significant & !is.na(result_df$log2FoldChange) & result_df$log2FoldChange > 0
  down <- significant & !is.na(result_df$log2FoldChange) & result_df$log2FoldChange < 0
  summary_rows[[i]] <- data.frame(
    contrast_id = contrast_id,
    label = label,
    factor = factor_name,
    numerator = numerator,
    denominator = denominator,
    tested_genes = sum(!is.na(result_df$pvalue)),
    significant = sum(significant),
    up = sum(up),
    down = sum(down),
    alpha = alpha
  )

  plot_df <- result_df
  plot_df$neg_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
  plot_df$status <- "Not significant"
  plot_df$status[up] <- "Up"
  plot_df$status[down] <- "Down"
  volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = status)) +
    geom_point(alpha = 0.55, size = 1.2, na.rm = TRUE) +
    scale_color_manual(values = c("Down" = "#2C7BB6", "Not significant" = "#BDBDBD", "Up" = "#D7191C")) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    theme_bw(base_size = 12) +
    labs(title = label, x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL)
  ggsave(file.path(contrast_dir, "volcano.png"), volcano, width = 8, height = 6, dpi = 180)

  ma <- ggplot(plot_df, aes(x = baseMean, y = log2FoldChange, color = status)) +
    geom_point(alpha = 0.55, size = 1.2, na.rm = TRUE) +
    scale_x_log10() +
    scale_color_manual(values = c("Down" = "#2C7BB6", "Not significant" = "#BDBDBD", "Up" = "#D7191C")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw(base_size = 12) +
    labs(title = paste(label, "MA plot"), x = "Mean normalized count", y = "log2 fold change", color = NULL)
  ggsave(file.path(contrast_dir, "ma.png"), ma, width = 8, height = 6, dpi = 180)
}

summary_table <- do.call(rbind, summary_rows)
write.table(
  summary_table,
  file.path(out_dir, "contrast_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("DESeq2 complete\n")
cat("Samples:", ncol(count_matrix), "\n")
cat("Tested genes:", nrow(count_matrix), "\n")
cat("Design:", deparse(design_formula), "\n")

