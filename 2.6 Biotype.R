library(dplyr)
library(edgeR)
library(ggplot2)
library(ggsci)
library(showtext)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

s <- catchSalmon(paths = file.path("/Users/MaxTomlinson/Desktop/Nanopore/Salmon", list.files(path = "/Users/MaxTomlinson/Desktop/Nanopore/Salmon")))

dge <- DGEList(counts = s$counts, genes = s$annotation)
dge.human <- dge[grep("^ENST", rownames(x = dge)), ]

vectors <- as.vector(x = rownames(x = dge.human$genes))

dge.human$genes$biotype <- sapply(X = strsplit(vectors, "\\|"), FUN = function(x) x[length(x = x)])

biotype_mapping <- c(
  "protein_coding_CDS_not_defined" = "protein_coding",
  "protein_coding_LoF" = "protein_coding",
  "pseudogene" = "pseudogene",
  "artifact" = "other",
  "non_stop_decay" = "other",
  "processed_transcript" = 'other',
  "ribozyme" = "other"
)

dge.human$genes$biotype <- sub(pattern = ".*pseudogene$", replacement = "pseudogene", x = dge.human$genes$biotype)
dge.human$genes$biotype <- sub(pattern = "^(IG|TR).*", replacement = "IG_or_TR_gene", x = dge.human$genes$biotype)
dge.human$genes$biotype <- sub(pattern = "nonsense_mediated_decay", replacement = "NMD", x = dge.human$genes$biotype)
dge.human$genes$biotype <- sub(pattern = "antisense", replacement = "lncRNA", x = dge.human$genes$biotype)
dge.human$genes$biotype <- sub(pattern = "non_coding", replacement = "lncRNA", x = dge.human$genes$biotype)
dge.human$genes$biotype <- sub(pattern = "retained_intron", replacement = "lncRNA", x = dge.human$genes$biotype)
dge.human$genes$biotype <- sub(pattern = "sense_intronic", replacement = "lncRNA", x = dge.human$genes$biotype)
dge.human$genes$biotype <- sub(pattern = "sense_overlapping", replacement = "lncRNA", x = dge.human$genes$biotype)
dge.human$genes$biotype <- ifelse(test = grepl(pattern = "RNA$", x = dge.human$genes$biotype) & !grepl(pattern = "lncRNA$", x = dge.human$genes$biotype), yes = "ncRNA", no = dge.human$genes$biotype)
dge.human$genes$biotype <- ifelse(test = dge.human$genes$biotype %in% names(x = biotype_mapping), yes = biotype_mapping[dge.human$genes$biotype], no = dge.human$genes$biotype)

data <- lapply(X = 1:18, FUN = function(x) {
  typesum <- aggregate(x = dge.human$counts[, x], by = list(dge.human$genes$biotype), FUN = sum, simplify = TRUE)
  return(typesum)
})

data <- do.call(what = "rbind", args = data)

sample <- sub(pattern = "_full_length_output", replacement = "", x = list.files(path = "/Users/MaxTomlinson/Desktop/Nanopore/Salmon"))

colnames(x = data) <- c("biotype", "total_count")

data$sample <- rep(x = sample, each = 8)
data$biotype <- gsub(pattern = "_", replacement = " ", x = data$biotype)
data$biotype <- gsub(pattern = "^(.)", replacement = "\\U\\1", x = data$biotype, perl = TRUE)
data$biotype <- gsub(pattern = "Nc", replacement = "nc", x = data$biotype)
data$biotype <- gsub(pattern = "Protein coding", replacement = "Protein-coding", x = data$biotype)

ord <- aggregate(x = data$total_count, by = list(data$biotype), FUN = sum, simplify = TRUE)
ord <- ord[order(ord$x, decreasing = TRUE), ]

plot <- ggplot(data = data, aes(x = sample, y = total_count, fill = factor(x = biotype, levels = ord$Group.1))) +
  geom_bar(stat = "identity", position = "fill") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Samples", y = "Proportion of transcript biotypes", fill = "Transcript biotype") +
  scale_fill_nejm() +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.6 Biotype.pdf", plot = plot, width = 11, height = 8.5)