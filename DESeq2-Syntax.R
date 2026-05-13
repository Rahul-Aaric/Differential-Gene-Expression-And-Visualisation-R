install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "EnhancedVolcano",
  "pheatmap"
))
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

counts <- read.csv("gene_counts.csv", row.names = 1)
metadata <- read.csv("metadata.csv", row.names = 1)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ batch + condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "case", "control"))

res <- lfcShrink(
  dds,
  contrast = c("condition", "case", "control"),
  res = res,
  type = "normal"
)

deg <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 1)

write.csv(as.data.frame(res), "DESeq2_all_results.csv")
write.csv(as.data.frame(deg), "DESeq2_significant_DEGs.csv")

#############################
# VISUALIZATION
#############################

pdf("MA_plot.pdf")
plotMA(res)
dev.off()

vsd <- vst(dds, blind = FALSE)

pdf("PCA_plot.pdf")
plotPCA(vsd, intgroup = "condition")
dev.off()

#############################
# VOLCANO PLOT
#############################

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

pdf("volcano_plot.pdf", width = 8, height = 6)

EnhancedVolcano(
  res_df,
  lab = res_df$gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Volcano Plot: Case vs Control",
  subtitle = "DESeq2 Differential Expression Results",
  xlab = "Log2 Fold Change",
  ylab = "-log10 adjusted p-value"
)

dev.off()
