library(dplyr)
library(ggplot2)
library(ggsci)
library(stringr)
library(tidyr)
library(tools)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

pull_metrics <- function(metric, lines) {
  line <- grep(pattern = metric, x = lines, value = TRUE)
  if (length(x = line) > 0) {
    parts <- str_split(string = line, pattern = "\\s+")[[1]]
    return(as.numeric(x = last(x = parts)))
  }
  return(NA)
}

process_file <- function(file_path) {
  lines <- readLines(con = file_path, warn = FALSE)
  
  tool <- file_path_sans_ext(x = basename(path = file_path))
  
  transcript_elements <- str_split(string = grep(pattern = "Transcript level:", x = lines, value = TRUE), pattern = "[[:space:]:|]+")[[1]]
  transcript_elements <- transcript_elements[transcript_elements != "" & !transcript_elements %in% c("Transcript", "level")]
  
  gene_elements <- str_split(string = grep(pattern = "Locus level:", x = lines, value = TRUE), pattern = "[[:space:]:|]+")[[1]]
  gene_elements <- gene_elements[gene_elements != "" & !gene_elements %in% c("Locus", "level")]
  
  sensitivity <- data.frame(
    Transcript = as.numeric(x = transcript_elements[1]),
    Gene = as.numeric(x = gene_elements[1]),
    Tool = tool
    )
  
  metrics <- c(
    "Matching Intron Chains" = "Matching intron chains:", 
    "Matching Transcripts" = "Matching transcripts:", 
    "Matching Loci" = "Matching loci:"
    )
  
  matching_metrics <- data.frame(
    Metric = names(x = metrics),
    Count = sapply(X = metrics, FUN = pull_metrics, lines),
    Tool = tool
    )
  
  features <- c(
    "Missed exons:", "Novel exons:", 
    "Missed introns:", "Novel introns:", 
    "Missed loci:", "Novel loci:"
    )
  
  feature_data <- lapply(X = features, FUN = function(pattern) {
    lines_matched <- grep(pattern = pattern, x = lines, value = TRUE)
    
    if (length(x = lines_matched) > 0) {
      feature <- str_remove(string = pattern, pattern = ":")
      type <- ifelse(test = str_detect(string = pattern, pattern = "Missed"), yes = "Missed", no = "Novel")
      counts <- as.numeric(x = str_extract(string = lines_matched, pattern = "\\d+/") %>% str_remove(pattern = "/"))
      percentage <- as.numeric(x = str_extract(string = lines_matched, pattern = "\\d+\\.?\\d*%") %>% str_remove(pattern = "%"))
      
      return(data.frame(Feature = feature, Type = type, Count = counts, Percentage = percentage, Tool = tool))
    } 
    
    return(data.frame(Feature = NA, Type = NA, Count = NA, Percentage = NA, Tool = tool))
  })
  
  feature_data <- do.call(what = rbind, args = feature_data) %>% arrange(Feature, Type)
  
  list(sensitivity = sensitivity, 
       matching_metrics = matching_metrics, 
       feature_data = feature_data)
}

file_paths <- list.files(path = "~/Desktop/Nanopore/Benchmarking/3.1+3.3 GffCompare", pattern = "\\.stats$", full.names = TRUE)

results <- lapply(X = file_paths, FUN = process_file)

all_sensitivity <- bind_rows(lapply(X = results, `[[`, "sensitivity")) %>%
  gather(key = "Metric", value = "Value", Transcript, Gene)
all_matching_metrics <- bind_rows(lapply(X = results, `[[`, "matching_metrics")) %>%
  group_by(Tool) %>%
  mutate(Percentage = Count / sum(Count) * 100)
all_feature_data <- bind_rows(lapply(X = results, `[[`, "feature_data"))

sensitivity_plot <- ggplot(data = all_sensitivity, aes(x = Tool, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(x = NULL, y = "Sensitivity (%)", fill = "Metric") +
  scale_fill_nejm(name = NULL) +
  scale_x_discrete(labels = c("Bambu", "FLAIR", "IsoQuant", "StringTie")) +
  theme(legend.position = "bottom")

matching_metrics_plot <- ggplot(data = all_matching_metrics, aes(x = Tool, y = Percentage, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(x = NULL, y = "Percentage", fill = "Metric") +
  scale_fill_nejm() +
  scale_x_discrete(labels = c("Bambu", "FLAIR", "IsoQuant", "StringTie"))

feature_data_plot <- ggplot(data = all_feature_data, aes(x = Tool, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Feature, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  scale_fill_nejm() +
  scale_x_discrete(labels = c("Bambu", "FLAIR", "IsoQuant", "StringTie")) +
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.1+3.3 GffCompare/Sensitivity.pdf", plot = sensitivity_plot, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.1+3.3 GffCompare/Matching metrics.pdf", plot = matching_metrics_plot, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.1+3.3 GffCompare/Feature data.pdf", plot = feature_data_plot, width = 11, height = 8.5)
