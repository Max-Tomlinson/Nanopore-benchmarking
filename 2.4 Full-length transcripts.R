library(dplyr)
library(GenomicFeatures)
library(ggplot2)
library(ggsci)
library(readxl)
library(rtracklayer)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none")
theme_set(new = theme_new)

RNA_extraction_data <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Nanopore/Extraction data/RNA_extraction_data_updated1[1].xlsx", skip = 1)

sample_info <- read.csv(file = "~/Desktop/Nanopore/sample_info.csv")
sample_info <- merge(x = sample_info, y = RNA_extraction_data[, c("...1", "RIN...13")], by.x = "No", by.y = "...1", sort = FALSE)

extract_info <- function(file_path) {
  lines <- readLines(con = file_path)
  sample_name <- gsub(pattern = ".fastq.gz", replacement = "", x = basename(path = lines[3]))
  percentage_of_reads <- as.numeric(x = regmatches(x = lines[length(x = lines) - 3], m = regexpr(pattern = "\\d+\\.\\d+", text = lines[length(x = lines) - 3])))
  data.frame(Sample = sample_name, Percentage = percentage_of_reads)
}

data <- do.call(what = rbind, lapply(X = list.files(path = "~/Desktop/Nanopore/Benchmarking/2.4 Full-length transcripts", pattern = "pychopper.out", full.names = TRUE, recursive = TRUE), FUN = extract_info))
data$Sample <- gsub(pattern = "sample_", replacement = "", x = data$Sample)
data <- cbind(data, RIN = sample_info$RIN...13[c(1:12, 14:19)])
data <- data[order(as.numeric(x = data$Sample)), ]

plot <- ggplot(data = data, aes(x = reorder(Sample, -Percentage), y = Percentage, fill = Percentage)) +
  coord_cartesian(ylim = c(75, 100)) +
  geom_bar(stat = "identity") + 
  labs(x = NULL, y = "Percentage of reads with two primers") +
  scale_fill_gradient(high = "#2999FF", low = "#FF0909", limits = c(min(data$Percentage), max(data$Percentage))) +
  theme(axis.text.x = element_blank())

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.4 Full-length transcripts/Barplot percentage of reads with two primers.pdf", plot = plot, width = 11, height = 8.5)

txdb <- makeTxDbFromGFF(file = "~/Desktop/Nanopore/References/gencode.v45.annotation.gtf", format = "gtf")

exons <- exonsBy(x = txdb, by = "tx", use.names = TRUE)

interval_length <- function(start, end) {
  return(end - start + 1)
}

data <- lapply(X = exons, FUN = function(exons_tx) {
  total_length <- sum(interval_length(start = start(x = exons_tx), end = end(x = exons_tx)))
  return(total_length)
})

data <- bind_rows(
  tx_id = names(x = data),
  tx_length = unlist(x = data)
)

density <- density(x = log(x = data$tx_length))

scaling <- max(density$y)

plot1 <- ggplot(data = data, aes(x = log(x = tx_length))) +
  geom_histogram(aes(y = after_stat(density) / scaling), bins = 100, fill = "#2999FF", color = "black", alpha = 0.5) +
  labs(x = "Log transcript length", y = "Proportion of full-length transcripts") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), sec.axis = sec_axis(~.*scaling, name = "Density"))

thresholds <- seq(min(log(x = data$tx_length)), max(log(x = data$tx_length)), length.out = 100)                                 

list <- list()

for (file_path in list.files(path = "~/Desktop/Nanopore/Benchmarking/2.4 Full-length transcripts", pattern = ".txt", full.names = TRUE, recursive = TRUE)) {
  data <- read.table(file = file_path, header = FALSE)
  data <- data.frame(Proportion = data$V1, Threshold = thresholds, Sample = gsub(pattern = "_full_length_output.sorted_coverage.txt", replacement = "", x = basename(path = file_path)))
  
  list[[file_path]] <- data
}

data <- do.call(what = rbind, args = list)

plot2 <- ggplot(data = data, aes(x = Threshold, y = Proportion, colour = factor(x = Sample)))
plot <- plot1 + 
  geom_line(data = ggplot_build(plot = plot2)$data[[1]], aes(x = x, y = y, colour = factor(x = data$Sample))) +
  guides(colour = "none")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.4 Full-length transcripts/Line plot percentage of full-length reads.pdf", plot = plot, width = 11, height = 8.5)

data <- merge(x = data, y = sample_info, by.x = "Sample", by.y = "Names", sort = FALSE)
data <- data[seq(100, nrow(x = data), by = 100), ]

cor_test <- cor.test(x = data$Proportion, y = data$RIN...13)

label_text <- paste("atop(italic(r) == ", round(x = cor_test$estimate, digits = 2), ", italic(p) == ", round(x = cor_test$p.value, digits = 2), ")")

plot <- ggplot(data = data, aes(x = Proportion, y = RIN...13, colour = factor(x = Sample))) +
  annotate("text", x = 0.70, y = 8.2, label = label_text, size = 9, parse = TRUE) +
  geom_smooth(method = "lm", se = TRUE, color = "#2999FF", fill = "#2999FF") +
  geom_point(size = 2, shape = 20, stroke = 5) +
  geom_point(colour = "black", size = 6, shape = 1, stroke = 1) +
  labs(x = "Proportion of full-length transcripts", y = "RIN")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.4 Full-length transcripts/Scatterplot full-length transcripts versus RIN.pdf", plot = plot, width = 11, height = 8.5)

BamSlam <- NULL

for (file in list.files(path = "~/Desktop/Nanopore/Benchmarking/2.4 Full-length transcripts/BamSlam", pattern = "\\.csv$", full.names = TRUE)) {
  BamSlam <- cbind(BamSlam, read.csv(file = file, header = TRUE)[, 2])
}

BamSlam <- t(x = BamSlam)

colnames(x = BamSlam) <- read.csv(file = file, header = TRUE)[, 1]

data <- cbind(data, BamSlam)

results <- data.frame()

for (col in colnames(x = data)[c(13:18, 21:ncol(x = data))]) {
  cor_test <- cor.test(x = data$RIN...13, y = data[[col]])
  
  results[col, "RIN_b"] <- round(x = cor_test$estimate, digits = 2)
  results[col, "RIN_p"] <- signif(x = cor_test$p.value, digits = 3)
  
  cor_test <- cor.test(x = data$Proportion, y = data[[col]])
  
  results[col, "Full_length_b"] <- round(x = cor_test$estimate, digits = 2)
  results[col, "Full_length_p"] <- signif(x = cor_test$p.value, digits = 3)
}

write.csv(x = results, "~/Desktop/Nanopore/Benchmarking/2.4 Full-length transcripts/Correlation.csv", row.names = TRUE)