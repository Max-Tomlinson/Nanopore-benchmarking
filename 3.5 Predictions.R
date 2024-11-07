library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(tidyr)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

gtf <- import(con = "~/Desktop/Nanopore/References/gencode.v45.annotation.gtf")

txdb <- makeTxDbFromGRanges(gr = gtf)

transcripts <- transcripts(x = txdb) %>%
  as.data.frame() %>%
  mutate(Annotation = "GENCODE", Type = "Transcript")

files <- list(
  Bambu = "~/Desktop/Nanopore/References/annotations_bambu.gtf",
  FLAIR = "~/Desktop/Nanopore/References/annotations_flair.gtf",
  IsoQuant = "~/Desktop/Nanopore/References/annotations_isoquant.gtf",
  StringTie = "~/Desktop/Nanopore/References/annotations_stringtie.gtf"
)

lengths <- data.frame()
for (tool in names(x = files)) {
  tool_txdb <- makeTxDbFromGFF(file = files[[tool]])
  tool_lengths <- transcripts(x = tool_txdb) %>%
    as.data.frame() %>%
    mutate(Annotation = tool, Type = "Transcript")
  lengths <- rbind(lengths, tool_lengths)
}

data <- merge(x = transcripts, y = lengths, by = c("tx_name", "Type"), all = TRUE)
data_plot <- data[complete.cases(data$width.y, data$width.x, data$Annotation.y), ]

scatter_plot <- ggplot(data = data_plot, aes(x = width.y, y = width.x, colour = Annotation.y)) +
  facet_wrap(~Annotation.y, scales = "free") +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed") +
  geom_point() +
  labs(x = "Expected transcript length", y = "Predicted transcript length") +
  scale_colour_nejm() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none")

gencode_ids <- unique(x = transcripts$tx_name)

data$Annotated <- ifelse(test = data$tx_name %in% gencode_ids, yes = "Yes", no = "No")

n_transcripts <- nrow(x = transcripts)

counts <- data %>%
  filter(Annotated == "Yes") %>%
  group_by(Annotation.y) %>%
  summarise(Count = n_distinct(tx_name), .groups = 'drop') %>%
  mutate(Proportion = Count / n_transcripts) %>%
  replace_na(replace = list(Annotation.y = "(Missing)"))

plot_data <- data %>%
  left_join(y = counts, by = c("Annotation.y")) %>%
  group_by(Annotation.y) %>%
  mutate(Total = n()) %>%
  group_by(Annotation.y, Annotated, Total) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  ungroup() %>%
  mutate(Proportion = Count / Total) %>%
  replace_na(list(Annotation.y = "(Missing)"))

bar_plot1 <- ggplot(data = plot_data, aes(x = Annotation.y, y = Count, fill = Annotated)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "Number of features") +
  scale_fill_nejm() +
  theme(legend.position = "bottom")

bar_plot2 <- ggplot(data = plot_data, aes(x = Annotation.y, y = Proportion, fill = Annotated)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "Proportion of features") +
  scale_fill_nejm() +
  theme(legend.position = "bottom")

bar_plot3 <- ggplot(data = counts[-5, ], aes(x = Annotation.y, y = Proportion, fill = Annotation.y)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = NULL, y = "Proportion of annotated transcripts found") +
  scale_fill_nejm() +
  theme(legend.position = "none")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.5 Predictions/Scatter.pdf", plot = scatter_plot, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.5 Predictions/Bar1.pdf", plot = bar_plot1, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.5 Predictions/Bar2.pdf", plot = bar_plot2, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.5 Predictions/Bar3.pdf", plot = bar_plot3, width = 11, height = 8.5)