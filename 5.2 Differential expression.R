library(DESeq2)
library(ggplot2)
library(gprofiler2)
library(ggsci)
library(IHW)
library(VennDiagram)

theme_new <-  theme_classic(base_size = 24) +
  theme(legend.position = "none", text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

colData <- read.csv(file = "~/Desktop/Nanopore/Analysis/Differential expression/colData.csv", row.names = 1)[1:18, ]

count_files <- c(
  "~/Desktop/Nanopore/References/bambu/counts_gene.txt", 
  "~/Desktop/Nanopore/References/FLAIR/gene_counts_ds2.tsv", 
  "~/Desktop/Nanopore/References/IsoQuant/OUT.gene_grouped_counts.tsv",
  "~/Desktop/Nanopore/References/stringtie/counts_gene.txt"
)
tools <- c("Bambu", "FLAIR", "IsoQuant", "StringTie")

DE_results <- data.frame(
  Tool = unique(x = tools),
  DE = rep(NA, length(x = unique(x = tools))),
  Up = rep(NA, length(x = unique(x = tools))),
  Down = rep(NA, length(x = unique(x = tools))))

data <- list()

for (i in seq_along(along.with = count_files)) {
  
  count_file <- count_files[i]
  
  tool <- tools[i]

  countData <- read.csv(file = count_file, sep = "\t", row.names = 1)
  
  colnames(x = countData) <- gsub(pattern = "^(.*?)_(.*?)_.*", replacement = "\\1_\\2", x = colnames(x = countData))
  
  if (tool == "FLAIR") {
    colnames(x = countData) <- gsub(pattern = "sample_([0-9])$", replacement = "sample_0\\1", x = gsub(pattern = "Sample", replacement = "sample_", x = gsub(pattern = "_.*", replacement = "", x = colnames(x = countData))))
  } else if (tool == "StringTie") {
    colnames(x = countData) <- gsub(pattern = "minimap2.", replacement = "", x = colnames(x = countData))
  }

  colData$Batch <- factor(x = colData$Batch)
  colData$Group <- factor(x = colData$Group)
  
  design <- as.formula("~ PC1 + PC2 + PC3 + RIN...13:Group + Batch + Group")
  
  dds <- DESeqDataSetFromMatrix(countData = round(x = countData), colData = colData, design = design)
  dds <- dds[rowSums(x = counts(object = dds) >= 10) >= 8, ]
  dds <- DESeq(object = dds)

  deseq2.results <- as.data.frame(x = DESeq2::results(object = dds, filterFun = ihw))
  deseq2.results <- deseq2.results[order(deseq2.results$padj), ]
  deseq2.results$Effect <- with(data = deseq2.results, ifelse(test = pvalue < 0.05 & log2FoldChange > 0, yes = "Higher", no = ifelse(test = pvalue < 0.05 & log2FoldChange < 0, yes = "Lower", no = NA)))
  deseq2.results$gene_id <- sub(pattern = "\\..*", replacement = "", x = rownames(x = deseq2.results))
  
  data[[tool]] <- deseq2.results
  
  DE <- sum(deseq2.results$pvalue <= 0.05, na.rm = TRUE)
  Up <- sum(deseq2.results$log2FoldChange > 0 & deseq2.results$pvalue <= 0.05, na.rm = TRUE)
  Down <- sum(deseq2.results$log2FoldChange < 0 & deseq2.results$pvalue <= 0.05, na.rm = TRUE)
  
  DE_results[DE_results$Tool == tool, "DE"] <- DE
  DE_results[DE_results$Tool == tool, "Up"] <- Up
  DE_results[DE_results$Tool == tool, "Down"] <- Down
  
}

data <- dplyr::bind_rows(data, .id = "tool")

write.csv(x = data, file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/DEGs.csv", row.names = FALSE)
write.csv(x = DE_results, file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/DEGs_summary.csv", row.names = FALSE)

plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  facet_wrap(~ tool, ncol = 2) +
  geom_point(data = data[is.na(data$Effect), ], aes(x = log2FoldChange, y = -log10(x = pvalue)), alpha = 0.25, size = 2/8, color = "#c1cad8") +
  geom_point(aes(colour = Effect), alpha = 0.5, size = 2, na.rm = TRUE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = expression(Log[2]~"fold change"), y = expression(-Log[10]~italic(P))) +
  scale_colour_nejm(name = "", breaks = c("Higher", "Lower")) +
  xlim(-max(abs(x = data$log2FoldChange)), max(abs(x = data$log2FoldChange)))

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Volcanos_gene.pdf", plot = plot, width = 11, height = 8.5)

gene_ids1 <- unique(x = data[which(x = data$tool == "Bambu" & data$pvalue < 0.05), ]$gene_id)
gene_ids2 <- unique(x = data[which(x = data$tool == "FLAIR" & data$pvalue < 0.05), ]$gene_id)
gene_ids3 <- unique(x = data[which(x = data$tool == "IsoQuant" & data$pvalue < 0.05), ]$gene_id)
gene_ids4 <- unique(x = data[which(x = data$tool == "StringTie" & data$pvalue < 0.05), ]$gene_id)

pdf(file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Venn_genes.pdf", width = 11, height = 8.5)

venn_plot <- venn.diagram(
  x = list(Bambu = gene_ids1, FLAIR = gene_ids2, IsoQuant = gene_ids3, StringTie = gene_ids4),
  category.names = c("Bambu", "FLAIR", "IsoQuant", "StringTie"),
  filename = NULL,
  output = TRUE,
  lwd = 2,
  lty = 'blank',
  col = c("#ba1f21", "#475da0", "#e27000", "#334429"),
  fill = c(scales::alpha("#ba1f21", 0.5), scales::alpha("#475da0", 0.5), scales::alpha("#e27000", 0.5), scales::alpha("#334429", 0.5)),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

grid.draw(x = venn_plot)

dev.off()
