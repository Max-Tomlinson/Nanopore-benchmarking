library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggsci)
library(IHW)
library(MatchIt)
library(patchwork)
library(scales)
library(tidyr)

theme_new <-  theme_classic(base_size = 24) + theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

pax_selected_samples <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/Nanopore/pax_selected_samples.txt", header = TRUE)
pax_selected_samples <- pax_selected_samples[!(pax_selected_samples$Study_No %in% c(371, 431, 1832, 2621, 32331, 371, 431, 32861, 51901)), 15]

colData1 <- read.csv(file = "~/Desktop/Nanopore/Analysis/Differential expression/colData.csv")[1:18, c(3:8)]
colData1 <- cbind(Study_No = 1:18, Age = pax_selected_samples, colData1)
colData1$Seq <- "Long"
colData2 <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/MetabPhenos_ALLEBSubjects_NoPhenoExclusions_WithHOMAs_20181120_dCell.txt", header = TRUE)[, c(1, 8, 69:71, 58, 63)]
colData2$Group <- with(data = colData2, ifelse(test = AGE_INCLUSION <= 55, yes = "Middle_aged", no = ifelse(test = AGE_INCLUSION >= 65, yes = "Older_aged", no = NA)))
colData2 <- colData2[which(x = colData2$BloodBatch >= 0), ]
colData2$Seq <- "Short"

colnames(x = colData2) <- colnames(x = colData1)

colData <- rbind(colData1, colData2)
colData <- na.omit(object = colData)
colData <- colData[order(colData$Age), ]
colData <- transform(`_data` = colData, Group = factor(x = Group), Batch = factor(x = Batch), Seq = factor(x = Seq))

middle_aged <- colData[which(x = colData$Group == "Middle_aged"), ]

matchit <- matchit(formula = Seq ~ Age + PC1 + PC2 + PC3, data = middle_aged, method = "cardinality")

middle_aged <- match.data(object = matchit, data = middle_aged)
middle_aged <- middle_aged[which(x = middle_aged$Seq == "Short"), ]

older_aged <- colData[which(x = colData$Group == "Older_aged"), ]

matchit <- matchit(formula = Seq ~ Age + PC1 + PC2 + PC3, data = older_aged, method = "cardinality")

older_aged <- match.data(object = matchit, data = older_aged)
older_aged <- older_aged[which(x = older_aged$Seq == "Short"), ]

colData <- rbind(middle_aged, older_aged)
colData[, 2:6] <- scale(x = colData[, 2:6])

countData <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Raw Counts/swap.merged.blood.gene.count.txt", header = TRUE, row.names = 4, check.names = FALSE)
countData <- countData[, colnames(x = countData) %in% colData$Study_No]

rownames(x = colData) <- colData$Study_No

design <- as.formula("~ PC1 + PC2 + PC3 + RIN...13 + Batch + Group")

dds <- DESeqDataSetFromMatrix(countData = round(x = countData), colData = colData, design = design)

keep <- rowSums(x = counts(object = dds) >= 10) >= 8

dds <- dds[keep,]
dds <- DESeq(object = dds)

deseq2.results <- as.data.frame(x = results(object = dds, filterFun = ihw))
deseq2.results$gene_id <- gsub(pattern = "\\..*", replacement = "", x = rownames(x = deseq2.results))

data <- read.csv(file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/DEGs.csv") %>%
  group_by(gene_id) %>%
  filter(n_distinct(tool) >= 4)

pvl_thresholds <- c(0.2, 0.1, 0.05, 0.01)
lfc_thresholds <- c(0, 0.5, 1, 1.5, 2)

results <- data.frame()

for (Tool in unique(x = data$tool)) {
  for (pvl in pvl_thresholds) {
    for (lfc in lfc_thresholds) {
      deseq2.results$Category <- ifelse(test = deseq2.results$pvalue < pvl & abs(x = deseq2.results$log2FoldChange) > lfc, yes = "positive", no = "negative")
      
      tool.results <- subset(x = data, tool == Tool)
      tool.results$Category <- ifelse(test = tool.results$pvalue < pvl & abs(x = tool.results$log2FoldChange) > lfc, yes = "positive", no = "negative")
      
      combined <- merge(x = deseq2.results, y = tool.results, by = "gene_id", suffixes = c(".short", ".long"))
      
      conf_matrix <- table(ShortRead = combined$Category.short, LongRead = combined$Category.long)
      
      Sensitivity <- round(x = conf_matrix["positive", "positive"] / sum(conf_matrix["positive", ]), digits = 3)
      Specificity <- round(x = conf_matrix["negative", "negative"] / sum(conf_matrix["negative", ]), digits = 3)
      
      results <- rbind(results, data.frame(Tool, pvl, lfc, Sensitivity, Specificity))
    }
  }
}

results <- pivot_wider(data = results, id_cols = c("pvl", "lfc"), names_from = "Tool", values_from = c("Sensitivity", "Specificity"), names_sep = "_")

write.csv(x = results, file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Short_read_sensitivity.csv", row.names = FALSE)

results <- pivot_longer(data = results, cols = -c(pvl, lfc), names_to = "variable", values_to = "value")
results <- results %>% separate(col = variable, into = c("metric", "tool"), sep = "_")

plot <- ggplot(results[which(x = results$pvl == 0.05), ], aes(x = lfc, y = value, color = tool, group = tool)) +
  geom_line(linewidth = 1) +
  facet_wrap(facets = ~ metric, ncol = 1, scales = "free_y") +
  labs(x = expression(Log[2]~"fold change threshold"), y = NULL) +
  scale_colour_nejm() + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_blank())

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Short_read_sensitivity.pdf", plot = plot, width = 11, height = 8.5, units = "in", useDingbats = FALSE)

results <- data.frame()

combined <- list()

for (Tool in unique(x = data$tool)) {
  tools <- subset(x = data, tool == Tool)
  
  combine <- merge(x = tools, y = deseq2.results, by = "gene_id", all = FALSE) %>%
    filter(abs(x = pvalue.y) < 1)
  overlap <- nrow(x = combine) / nrow(x = tools) * 100
  
  same_effect <- sum(sign(x = combine$log2FoldChange.x) == sign(x = combine$log2FoldChange.y), na.rm = TRUE)
  rate_effect <- same_effect / nrow(x = combine) * 100
  both_significant <- sum(combine$pvalue.x < 0.05 & combine$pvalue.y < 0.05, na.rm = TRUE)
  rate_significant <- both_significant / sum(combine$pvalue.y < 0.05) * 100
  
  cor_test <- cor.test(x = abs(x = combine[combine$pvalue.y < 0.05, ]$log2FoldChange.x), y = abs(x = combine[combine$pvalue.y < 0.05, ]$log2FoldChange.y))
  
  results <- rbind(results, data.frame(
    Tool            = Tool,
    Overlap         = nrow(x = combine),
    RateOverlap     = signif(x = overlap, digits = 3),
    SameEffect      = same_effect,
    RateEffect      = signif(x = rate_effect, digits = 3),
    BothSignificant = both_significant,
    RateSignificant = signif(x = rate_significant, digits = 3),
    CorEstimate     = signif(x = cor_test$estimate, digits = 3),
    CorPvalue       = signif(x = cor_test$p.value, digits = 3)
  ))
  
  combined[[Tool]] <- combine
}

write.csv(x = results, file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Short_read_overlap.csv", row.names = FALSE)

data <- pivot_longer(data = results, cols = c("RateOverlap", "RateEffect", "RateSignificant"), names_to = "Result", values_to = "Rate")
data$Result <- factor(x = data$Result, levels = c("RateOverlap", "RateEffect", "RateSignificant"))

plot1 <- ggplot(data = data, aes(x = Tool, y = Rate, fill = Result)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") +
  labs(fill = NULL, x = NULL, y = "Overlap rate with short-read data") +
  scale_fill_nejm(labels = c("Overlap", "Effect", "Significant")) + 
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(size = 16), legend.position = "bottom", text = element_text(family = "Helvetica"))

data <- do.call(what = rbind, args = combined)

plot2 <- ggplot(data[data$pvalue.y < 0.05, ], aes(x = abs(x = log2FoldChange.x), y = abs(x = log2FoldChange.y), colour = tool, fill = tool)) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(colour = NULL, fill = NULL, x = expression(Log[2]~"fold change (long)"), y = expression(Log[2]~"fold change (short)")) +
  scale_colour_nejm() + 
  scale_fill_nejm() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  theme_classic(base_size = 24) +
  theme(legend.position = "bottom", legend.text = element_text(size = 11), text = element_text(family = "Helvetica"))

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Short_read_overlap.pdf", plot = plot1, width = 11, height = 8.5, units = "in", useDingbats = FALSE)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Short_read_correlation.pdf", plot = plot2, width = 11, height = 8.5, units = "in", useDingbats = FALSE)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Short_read_combined.pdf", plot = (plot + plot2), width = 11, height = 8.5, units = "in", useDingbats = FALSE)