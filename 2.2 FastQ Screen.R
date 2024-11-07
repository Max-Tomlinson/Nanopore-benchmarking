library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(readr)
library(tidyr)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

data_files <- list.files(path = "~/Desktop/Nanopore/Benchmarking/2.2 FastQ Screen", pattern = ".txt", full.names = TRUE, recursive = TRUE)
genomes <- c("Human", "Mouse", "Rat", "Drosophila", "Worm", "Yeast", "Ecoli", "rRNA", "MT", "Vectors")

process_genome <- function(genome) {
  genome_data <- bind_rows(lapply(X = data_files, FUN = function(file) {
    sample_data <- read_delim(file = file, col_names = TRUE, skip = 1, delim = "\t", n_max = 14, trim_ws = TRUE, show_col_types = FALSE)
    sample_data$Sample <- basename(path = file)
    sample_data
  }))
  
  genome_subset   <- filter(.data = genome_data, Genome == genome)
  target_columns  <- c("%One_hit_one_genome", "%One_hit_multiple_genomes", "%Multiple_hits_one_genome", "%Multiple_hits_multiple_genomes")
  Mean_Percentage <- colMeans(x = genome_subset[, target_columns], na.rm = TRUE)
  
  genome_summary <- data.frame(Genome = genome, target_columns, Mean_Percentage = Mean_Percentage)
  return(genome_summary)
}

data <- lapply(X = genomes, FUN = process_genome)
data <- do.call(what = rbind, args = data)
data$Alignment_Status <- gsub(pattern = "%", replacement = "", x = data$target_columns)
data$Alignment_Status <- gsub(pattern = "_", replacement = " ", x = data$Alignment_Status)
data$Genome <- factor(x = data$Genome, levels = genomes, ordered = TRUE)

plot <- ggplot(data = data, aes(x = Genome, y = Mean_Percentage, fill = Alignment_Status)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  labs(fill = "Alignment status", x = NULL, y = "Average percentage identity") +
  scale_fill_nejm() +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 16))
  
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.2 FastQ Screen/FastQ Screen plot.pdf", plot = plot, width = 11, height = 8.5, units = "in", useDingbats = FALSE)
