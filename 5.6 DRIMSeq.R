library(DRIMSeq)
library(ggplot2)
library(gprofiler2)
library(stageR)

theme_new <- theme_classic(base_size = 24) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

samples <- read.csv(file = "~/Desktop/Nanopore/Analysis/Differential expression/colData.csv")[1:18, ]
samples$Batch <- factor(x = samples$Batch)
samples$Group <- factor(x = samples$Group)

colnames(x = samples)[c(1, 8)] <- c("sample_id", "condition")

counts <- data.frame(read.table(file = "~/Desktop/Nanopore/References/bambu/counts_transcript.txt", header = TRUE, sep = "\t"))

colnames(x = counts) <- sub(pattern = "^(.*?)_(.*?)_.*", "\\1_\\2", x = colnames(x = counts))
colnames(x = counts)[1:2] <- c("feature_id", "gene_id")

d <- dmDSdata(counts = counts, samples = samples)
d <- dmFilter(x = d, min_samps_gene_expr = 18, min_samps_feature_expr = 8, min_gene_expr = 10, min_feature_expr = 10)

design <- model.matrix(~ Age + PC1 + PC2 + PC3 + RIN...13:condition + Batch + condition, data = DRIMSeq::samples(x = d))

set.seed(seed = 123)

d <- dmPrecision(x = d, design = design)
d <- dmFit(x = d, design = design)
d <- dmTest(x = d, coef = "conditionOlder_aged")

res <- DRIMSeq::results(x = d)
res <- res[order(res$pvalue), ]

Precision   <- plotPrecision(x = d) + labs(x = bquote(Log[10] ~ "of mean expression"), y = bquote(Log[10] ~ "precision")) + theme_classic(base_size = 24) + theme(legend.position = "top", text = element_text(family = "Helvetica"))
PValues     <- plotPValues(x = d, level = "gene") + labs(x = "P-values") + theme_classic(base_size = 24) + theme(text = element_text(family = "Helvetica"))
Proportions <- plotProportions(x = d, gene_id = res$gene_id[1], group_variable = "condition", plot_type = "boxplot1") + labs(title = "TSEN34", x = NULL) + theme_new + theme(axis.text.x = element_text(size = 14))

pScreen <- DRIMSeq::results(x = d)$pvalue
names(x = pScreen) <- DRIMSeq::results(x = d)$gene_id
pScreen <- na.omit(object = pScreen)

pConfirmation <- matrix(data = DRIMSeq::results(x = d, level = "feature")$pvalue, ncol = 1)
rownames(x = pConfirmation) <- DRIMSeq::results(x = d, level = "feature")$feature_id
pConfirmation <- na.omit(object = pConfirmation)

tx2gene <- DRIMSeq::results(x = d, level = "feature")[, c("feature_id", "gene_id")]

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu", alpha = 0.1, allowNA = FALSE)

SignificantGenes <- getSignificantGenes(object = stageRObj)
SignificantTX    <- getSignificantTx(object = stageRObj)

padj <- getAdjustedPValues(object = stageRObj, order = TRUE, onlySignificantGenes = FALSE)

pConfirmation <- pConfirmation[order(pConfirmation), ]

query <- unique(x = sub(pattern = "\\..*", replacement = "", x = names(x = pConfirmation[startsWith(x = names(x = pConfirmation), prefix = "ENS") & pConfirmation < 0.05])))

gostres <- gost(query = query, organism = "hsapiens", ordered_query = TRUE, significant = TRUE)

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Precision.pdf", plot = Precision, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/PValues.pdf", plot = PValues, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Proportions.pdf", plot = Proportions, width = 11, height = 8.5)

write.csv(x = as.data.frame(x = gostres$result)[, -14], file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/gostres_DRIMSeq.csv", row.names = FALSE)
write.csv(x = res, file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/DRIMSeq.csv", row.names = FALSE)