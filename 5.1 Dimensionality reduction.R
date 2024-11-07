library(DESeq2)
library(ggplot2)
library(ggsci)
library(rrcov)
library(umap)

theme_new <-  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

count_files <- c(
  "~/Desktop/Nanopore/References/bambu/counts_gene.txt", 
  "~/Desktop/Nanopore/References/FLAIR/gene_counts_ds2.tsv", 
  "~/Desktop/Nanopore/References/IsoQuant/OUT.gene_grouped_counts.tsv",
  "~/Desktop/Nanopore/References/stringtie/counts_gene.txt"
)

tools <- c("Bambu", "FLAIR", "IsoQuant", "StringTie")

plot_data_pca <- list()
plot_data_umap <- list()
plot_data_rpca <- list()

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
  
  colData <- data.frame(colnames = colnames(x = countData))
  colData$condition <- factor(x = ifelse(test = grepl(pattern = "^sample_0[4-9]|^sample_1[0-1]", x = colData$colnames), yes = "Middle_aged", no = "Older_aged"))
  
  dds <- DESeqDataSetFromMatrix(countData = round(x = countData), colData = colData, design = as.formula(object = "~ condition"))
  
  keep <- rowSums(x = counts(object = dds) >= 10) >= 8
  
  dds <- dds[keep,]
  dds <- estimateSizeFactors(object = dds)
  
  rld <- rlog(object = dds, blind = FALSE)
  rld$condition <- gsub(pattern = "_", replacement = "-", x = rld$condition)
  
  pcaData <- plotPCA(object = rld, intgroup = c("condition"), returnData = TRUE)
  pcaData$Tool <- tool
  
  umapData <- umap(d = t(x = assay(x = rld)))
  umapData <- data.frame(UMAP1 = umapData$layout[, 1], UMAP2 = umapData$layout[, 2], condition = colData(rld)$condition, Tool = tool)
  umapData$name <- rownames(x = umapData)
  
  robustPCAData <- PcaGrid(x = as.matrix(x = t(x = assay(x = rld))))
  robustPCAData <- data.frame(`Score distance` = robustPCAData$sd, `Score distance cut-off` = robustPCAData$cutoff.sd, Index = colnames(x = rld), condition = colData(rld)$condition, Tool = tool)
  
  plot_data_pca[[length(plot_data_pca) + 1]] <- pcaData
  plot_data_umap[[length(plot_data_umap) + 1]] <- umapData
  plot_data_rpca[[length(plot_data_rpca) + 1]] <- robustPCAData
}

plot_data <- do.call(what = rbind, args = plot_data_pca)

pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, colour = condition, shape = name)) +
  facet_wrap(~ Tool, nrow = 2, ncol = 2) +
  geom_point(size = 3) +
  guides(colour = guide_legend(title = NULL), shape = "none") +
  scale_colour_nejm() +
  scale_shape_manual(values = 1:19)

plot_data <- do.call(what = rbind, args = plot_data_umap)

umap_plot <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, colour = condition, shape = name)) +
  facet_wrap(~ Tool, nrow = 2, ncol = 2) +
  geom_point(size = 3) +
  guides(colour = guide_legend(title = NULL), shape = "none") +
  scale_colour_nejm() +
  scale_shape_manual(values = 1:19)

plot_data <- do.call(what = rbind, args = plot_data_rpca)

rpc_plot <- ggplot(plot_data, aes(x = factor(x = Index), y = Score.distance, colour = condition, shape = Index)) +
  facet_wrap(~ Tool, nrow = 2, ncol = 2) +
  geom_hline(yintercept = mean(x = plot_data$Score.distance.cut.off), colour = "red") +
  geom_point(size = 3) +
  guides(colour = guide_legend(title = NULL), shape = "none") +
  scale_colour_nejm() +
  scale_shape_manual(values = 1:19) +
  theme(axis.text.x = element_blank(), axis.title = element_blank())

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.1 Dimensionality reduction/PCA.pdf", plot = pca_plot, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.1 Dimensionality reduction/UMAP.pdf", plot = umap_plot, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.1 Dimensionality reduction/RPC.pdf", plot = rpc_plot, width = 11, height = 8.5)