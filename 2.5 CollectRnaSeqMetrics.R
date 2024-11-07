library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ggsci)
library(scales)
library(tidyr)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

sample_info <- read.csv(file = "~/Desktop/Nanopore/sample_info.csv", row.names = 9)

coverages <- NULL
alignment <- NULL

for (file_path in list.files(path = "~/Desktop/Nanopore/Benchmarking/2.5 CollectRnaSeqMetrics", pattern = "Metrics", full.names = TRUE)) {
  metrics <- read.table(file = file_path, header = TRUE, skip = 8)
  histogram <- read.table(file = file_path, header = FALSE, skip = 7, nrows = 1)
  coverages <- cbind(coverages, metrics$All_Reads.normalized_coverage)
  alignment <- rbind(alignment, histogram)
}

data <- as.data.frame(x = coverages)
colnames(x = data) <- list.files(path = "~/Desktop/Nanopore/Benchmarking/2.5 CollectRnaSeqMetrics", pattern = "Metrics", all.files = TRUE, recursive = TRUE)
colnames(x = data) <- gsub(pattern = ".metrics", replacement = "", x = colnames(x = data))

data <- gather(data = data, key = "Sample", value = "NormalisedCoverage", factor_key = TRUE)
data$`Normalised position` <- as.numeric(x = rep(x = metrics$normalized_position, 18))

plot <- ggplot(data = data, aes(x = `Normalised position`, y = `NormalisedCoverage`, colour = Sample)) +
  geom_line() +
  guides(colour = "none") +
  labs(x = "Normalised position", y = "Normalised coverage") +
  scale_fill_nejm()

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.5 CollectRnaSeqMetrics/Normalised gene coverage.pdf", plot = plot, width = 11, height = 8.5)

data <- as.data.frame(x = alignment, row.names = list.files(path = "~/Desktop/Nanopore/Benchmarking/2.5 CollectRnaSeqMetrics", pattern = "Metrics", all.files = TRUE, full.names = FALSE, recursive = TRUE))
colnames(x = data) <- read.table(file = file_path, header = FALSE, skip = 6, nrows = 1)[, -c(3, 16, 28:30)]

data <- gather(data = data[, c(3:6)], key = "Metric", value = "Value", factor_key = TRUE)
replacements <- c("CODING_BASES" = "Coding bases", "UTR_BASES" = "UTR bases", "INTRONIC_BASES" = "Intronic bases", "INTERGENIC_BASES" = "Intergenic bases")
data$Metric <- replacements[data$Metric]
data$Sample <- list.files(path = "~/Desktop/Nanopore/Benchmarking/2.5 CollectRnaSeqMetrics", pattern = "Metrics")

plot <- ggplot(data = data, aes(x = Sample, y = Value, fill = Metric)) + 
  coord_flip() +
  geom_bar(stat = "identity", position = position_fill(), width = 0.5) +
  labs(y = "Proportion") + 
  scale_fill_nejm() +
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 18))

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.5 CollectRnaSeqMetrics/Alignment summary.pdf", plot = plot, width = 11, height = 8.5)