library(DEXSeq)

dxr <- readRDS(file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/dxr.rds")

pdf(file = "~/Desktop/Nanopore/Benchmarking/5.2 Differential Expression/DEXSeq_LRRK2.pdf", width = 11.5, height = 8)

plotDEXSeq(object = dxr, geneID = "ENSG00000188906.17", expression = FALSE, splicing = TRUE, displayTranscripts = FALSE, legend = TRUE, color = c("#4040FF", "#FF4162"), cex.axis = 1.2, cex = 2, lwd = 3)

dev.off()
