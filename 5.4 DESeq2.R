library(EnhancedVolcano)
library(DESeq2)
library(gprofiler2)
library(IHW)

colData <- read.csv(file = "~/Desktop/Nanopore/Analysis/Differential expression/colData.csv", row.names = 1)[1:18, ]
colData$Batch <- factor(x = colData$Batch)
colData$Group <- factor(x = colData$Group)

design <- as.formula("~ PC1 + PC2 + PC3 + RIN...13:Group + Batch + Group")

countData <- read.csv(file = "~/Desktop/Nanopore/References/bambu/counts_transcript.txt", sep = "\t", row.names = 1)[, -1]

colnames(x = countData) <- sub(pattern = "^(.*?)_(.*?)_.*", replacement = "\\1_\\2", x = colnames(x = countData))

dds <- DESeqDataSetFromMatrix(countData = round(x = countData), colData = colData, design = design)

keep <- rowSums(x = counts(object = dds) >= 10) >= 8

dds <- dds[keep,]
dds <- DESeq(object = dds)

res <- DESeq2::results(object = dds, filterFun = ihw)
res <- as.data.frame(x = res[order(res$padj), ])
res$transcript_id <- sub(pattern = "\\..*", replacement = "", x = rownames(x = res))

transcript_ids <- gconvert(query = res$transcript_id, organism = "hsapiens", mthreshold = 1, filter_na = TRUE)

res$names <- ifelse(test = !is.na(x = match(x = res$transcript_id, table = transcript_ids$input)), yes = transcript_ids$name[match(x = res$transcript_id, table = transcript_ids$input)], no = res$transcript_id)
res$names[7] <- "AC087392.4"

gostres <- gost(query = res[which(x = startsWith(x = res$transcript_id, prefix = "ENS") & res$pvalue < 0.05), ]$transcript_id, organism = "hsapiens", ordered_query = TRUE)

write.csv(x = as.data.frame(x = gostres$result)[, -14], file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/gostres_bambu.csv", row.names = FALSE)
write.csv(x = res, file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/DESeq.csv", row.names = FALSE)

lab_italics <- paste0("italic('", make.unique(names = res$names), "')")

keyvals <- ifelse(test = res$log2FoldChange < 0 & res$padj < 0.1, yes = '#475da0', no = ifelse(test = res$log2FoldChange > 0 & res$padj < 0.1, yes = '#ba1f21', no = 'grey'))

names(x = keyvals)[keyvals == '#ba1f21'] <- 'Up-regulated'
names(x = keyvals)[keyvals == 'grey'] <- 'NS'
names(x = keyvals)[keyvals == '#475da0'] <- 'Down-regulated'

plot <- EnhancedVolcano(toptable = res,
                        lab = lab_italics,
                        x = 'log2FoldChange',
                        y = 'padj',
                        selectLab = lab_italics[1:20],
                        axisLabSize = 24,
                        title = NULL,
                        subtitle = NULL,
                        caption = NULL,
                        pCutoff = 0.1,
                        pCutoffCol = 'padj',
                        FCcutoff = 0,
                        pointSize = c(ifelse(test = abs(x = res$log2FoldChange) > 0 & res$padj < 0.1, yes = 2, no = 2/8)),
                        boxedLabels = TRUE,
                        parseLabels = TRUE,
                        colCustom = keyvals,
                        legendPosition = "none",
                        legendLabSize = 24,
                        drawConnectors = TRUE,
                        arrowheads = FALSE,
                        max.overlaps = 30,
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE)

ggsave("~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/Volcano_bambu.pdf", plot = plot, width = 11, height = 8.5, useDingbats = FALSE)
