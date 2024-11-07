library(dplyr)
library(ggplot2)
library(ggsci)
library(readxl)
library(tidyr)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

read_custom_data <- function(file_path, num_lines) {
  lines <- readLines(con = file_path, n = num_lines)
  processed_lines <- lapply(X = lines, FUN = function(line) {
    parts <- strsplit(x = line, split = "\t", fixed = TRUE)[[1]]
  })
  data_frame <- do.call(what = rbind, args = processed_lines)
  data_frame <- as.data.frame(x = data_frame)
  colnames(x = data_frame) <- data_frame[1, ]
  data_frame <- data_frame[-1, ]
  data_frame <- data_frame[, order(colnames(x = data_frame))]
  return(data_frame)
}

fastq_files <- read_custom_data(file_path = file.path("~/Desktop/Nanopore/Benchmarking/2.3 NanoComp/FastQ_NanoStats.txt"), num_lines = 9)
fastq_files <- pivot_longer(data = fastq_files, cols = -Metrics, names_to = "Sample", values_to = "FASTQ")  %>% mutate(FASTQ = as.numeric(FASTQ))

bam_files <- read_custom_data(file_path = file.path("~/Desktop/Nanopore/Benchmarking/2.3 NanoComp/BAM_NanoStats.txt"), num_lines = 13)
bam_files <- pivot_longer(data = bam_files, cols = -Metrics, names_to = "Sample", values_to = "BAM")  %>% mutate(BAM = as.numeric(BAM))

data <- merge(x = fastq_files, y = bam_files, by = c("Metrics", "Sample"))
data$Metrics <- sub(pattern = "^(.)", replacement = "\\U\\1", x = gsub(pattern = "qual", replacement = "quality", x = gsub(pattern = "n50", replacement = "N50", x = gsub(pattern = "stdev", replacement = "STDEV", x = gsub(pattern = "_", replacement = " ", x = data$Metrics)))), perl = TRUE)

plots <- lapply(X = split(x = data, f = data$Metric), FUN = function(data) {
  plot <- ggplot(data = data, aes(x = FASTQ, y = BAM, fill = factor(x = Sample))) +
    geom_point(size = 5, shape = 21, color = "black", alpha = 0.8) +
    geom_abline(intercept = 0, color = "grey70", linetype = "dashed") +
    labs(x = paste(data$Metrics[1], "(FastQ)"), y = paste(data$Metrics[1], "(BAM)")) +
    theme(legend.position = "none")
  plot
})

for (i in names(x = plots)) {
  file_name <- paste0("~/Desktop/Nanopore/Benchmarking/2.3 NanoComp/", i, ".pdf")
  ggsave(filename = file_name, plot = plots[[i]], width = 11, height = 8.5)
}

RNA_extraction_data <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Nanopore/Extraction data/RNA_extraction_data_updated1[1].xlsx", skip = 1)

sample_info <- read.csv(file = "~/Desktop/Nanopore/sample_info.csv")[c(1:12, 14:19), ]
sample_info <- merge(x = sample_info, y = RNA_extraction_data[, c("...1", "RIN...13")], by.x = "No", by.y = "...1", sort = FALSE)

data <- merge(x = data, y = sample_info[, c("Names", "RIN...13")], by.x = "Sample", by.y = "Names", sort = FALSE)

perform_cor_test <- function(metric_data) {
  cor_test_raw <- cor.test(x = metric_data$FASTQ, y = metric_data$RIN...13, method = "pearson")
  cor_test_bam <- cor.test(x = metric_data$BAM, y = metric_data$RIN...13, method = "pearson")
  
  list(mean_raw = mean(metric_data$FASTQ),
       cor_raw = cor_test_raw$estimate,
       p_value_raw = cor_test_raw$p.value,
       mean_bam = mean(metric_data$BAM),
       cor_bam = cor_test_bam$estimate,
       p_value_bam = cor_test_bam$p.value)
}

results <- data %>%
  group_by(Metrics) %>%
  do(data.frame(perform_cor_test(.))) %>%
  mutate(across(where(is.numeric), ~signif(., digits = 3)))

write.csv(x = results, "~/Desktop/Nanopore/Benchmarking/2.3 NanoComp/Correlation_RIN.csv", row.names = FALSE)

average_read_lengths <- read.table(file = "~/Desktop/Nanopore/Benchmarking/2.3 NanoComp/average_read_lengths.txt", check.names = FALSE)
average_read_lengths <- average_read_lengths[order(average_read_lengths$V1), ]

data <- data[data$Metric == "Mean read length", ]
data <- cbind(data, GENCODE = average_read_lengths$V2)

plot <- ggplot(data = data, aes(x = BAM, y = GENCODE, fill = factor(x = Sample))) +
  geom_point(size = 5, shape = 21, color = "black", alpha = 0.8) +
  geom_abline(intercept = 0, color = "grey70", linetype = "dashed") +
  labs(x = "Mean read length (BAM)", y = "Mean read length (GENCODE)") +
  theme(legend.position = "none")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.3 NanoComp/Mean read length (GENCODE).pdf", plot = plot, width = 11, height = 8.5)