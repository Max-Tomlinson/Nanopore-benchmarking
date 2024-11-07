library(BSgenome.Hsapiens.NCBI.GRCh38)
library(edgeR)
library(GenomicFeatures)
library(IsoformSwitchAnalyzeR)
library(tximeta)

isoformCountMatrix <- read.table(file = "~/Desktop/Nanopore/References/bambu/counts_transcript.txt", header = TRUE, sep = "\t", row.names = 1)[, -1]

myDesign <- data.frame(
  sampleID  = colnames(x = isoformCountMatrix),
  condition = factor(x = ifelse(test = grepl("^sample_0[4-9]|^sample_1[0-1]", x = colnames(x = isoformCountMatrix)), yes = "Middle_aged", no = "Older_aged")),
  batch     = c(3, 3, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
)

mySwitchList <- importRdata(
  isoformCountMatrix            = isoformCountMatrix,
  isoformExonAnnoation          = "~/Desktop/Nanopore/References/bambu/extended_annotations.gtf",
  designMatrix                  = myDesign,
  addAnnotatedORFs              = FALSE,
  removeTECgenes                = FALSE,
  estimateDifferentialGeneRange = FALSE
)

mySwitchList <- preFilter(switchAnalyzeRlist = mySwitchList, geneExpressionCutoff = NULL, IFcutoff = NULL)
mySwitchList <- isoformSwitchTestSatuRn(switchAnalyzeRlist = mySwitchList, reduceToSwitchingGenes = FALSE)
mySwitchList <- analyzeAlternativeSplicing(switchAnalyzeRlist = mySwitchList, alpha = 0.1)

localTheme <- theme_classic() + theme(text = element_text(family = "Helvetica"))

Summary <- extractSplicingSummary(switchAnalyzeRlist = mySwitchList, localTheme = localTheme)
Enrichment <- extractSplicingEnrichment(switchAnalyzeRlist = mySwitchList, localTheme = localTheme, returnResult = TRUE, minEventsForPlotting = 1)
GenomeWide <- extractSplicingGenomeWide(switchAnalyzeRlist = mySwitchList, localTheme = localTheme, returnResult = FALSE)

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/IsoformSwitchAnalyzeR/Summary.pdf", plot = Summary, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/IsoformSwitchAnalyzeR/GenomeWide.pdf", plot = GenomeWide, width = 11, height = 8.5)

Enrichment$AStype <- c("A3", "A5", "AL", "AF", "SE", "IR", "MX", "MSE")
Enrichment$size <- with(data = Enrichment, expr = nUp + nDown)

plot <- ggplot(Enrichment, aes(x = propUp, y = reorder(AStype, -propUp))) +
  geom_errorbarh(aes(xmin = propUpCiLo, xmax = propUpCiHi, colour = Significant), height = 0.2) +
  geom_point(aes(colour = Significant, size = size)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
  labs(x = "Enrichment of alternative splicing events", y = NULL, colour = "FDR 5%", size = "Genes") +
  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica")) +
  scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'black'), labels = c("False", "True")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/IsoformSwitchAnalyzeR/Enrichment.pdf", plot = plot, width = 11, height = 8.5)

# stringTieQuant <- importIsoformExpression(
#   parentDir  = "~/Desktop/Nanopore/Ballgown",
#   readLength = 583,
# )
#
# myDesign <- data.frame(
#   sampleID  = colnames(x = stringTieQuant$abundance)[2:19],
#   condition = factor(x = ifelse(test = grepl("^sample_0[4-9]|^sample_1[0-1]", x = colnames(x = stringTieQuant$abundance)[2:19]), yes = "Middle_aged", no = "Older_aged")),
#   batch     = c(3, 3, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
# )
# 
# mySwitchList <- importRdata(
#   isoformCountMatrix            = stringTieQuant$counts,
#   isoformRepExpression          = stringTieQuant$abundance,
#   designMatrix                  = myDesign,
#   isoformExonAnnoation          = "~/Desktop/Nanopore/References/stringtie/stringtie.merged.gtf",
#   detectUnwantedEffects         = FALSE,
#   addAnnotatedORFs              = FALSE,
#   removeTECgenes                = FALSE,
#   estimateDifferentialGeneRange = FALSE
# )
# 
# salmonQuant <- importIsoformExpression(
#   parentDir  = "~/Desktop/Nanopore/Salmon"
# )
# 
# myDesign <- data.frame(
#   sampleID  = colnames(x = salmonQuant$abundance)[2:19],
#   condition = factor(x = ifelse(test = grepl("^sample_0[4-9]|^sample_1[0-1]", x = colnames(x = salmonQuant$abundance)[2:19]), yes = "Middle_aged", no = "Older_aged")),
#   batch     = c(3, 3, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
# )
# 
# mySwitchList <- importRdata(
#   isoformCountMatrix            = salmonQuant$counts,
#   isoformRepExpression          = salmonQuant$abundance,
#   isoformExonAnnoation          = "~/Desktop/Nanopore/References/gencode.v45.annotation.gtf",
#   isoformNtFasta                = "~/Desktop/Nanopore/References/gencode.v45.transcripts.fa",
#   designMatrix                  = myDesign,
#   detectUnwantedEffects         = FALSE,
#   addAnnotatedORFs              = FALSE,
#   removeTECgenes                = FALSE,
#   estimateDifferentialGeneRange = FALSE
# )
