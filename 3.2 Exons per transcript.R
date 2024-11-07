library(dplyr)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(RColorBrewer)
library(rtracklayer)
library(tidyr)

gtf_files <- c(
  "~/Desktop/Nanopore/References/annotations_bambu.gtf",
  "~/Desktop/Nanopore/References/annotations_flair.gtf",
  "~/Desktop/Nanopore/References/annotations_isoquant.gtf",
  "~/Desktop/Nanopore/References/annotations_stringtie.gtf")

tools <- c("Bambu", "FLAIR", "IsoQuant", "StringTie")

all_data_TX <- data.frame(Transcript = character(), Exons = integer(), Tool = character())
all_data_GN <- data.frame(Gene = character(), Transcripts = integer(), Tool = character())

for (i in seq_along(along.with = gtf_files)) {

  gtf <- import(con = gtf_files[i])
  
  exons <- gtf[gtf$type == "exon"]
  exons_per_transcript <- table(exons$transcript_id)
  exons_per_transcript <- as.data.frame(x = exons_per_transcript, responseName = "Exons")
  
  names(x = exons_per_transcript) <- c("Transcript", "Exons")
  
  exons_per_transcript$Tool <- tools[i]
  
  transcript <- gtf[gtf$type == "transcript"]
  
  if (i == 4){
    transcripts_per_gene <- table(transcript$ref_gene_id)
  } else {
    transcripts_per_gene <- table(transcript$gene_id)
  }
  
  transcripts_per_gene <- as.data.frame(x = transcripts_per_gene, responseName = "Transcripts")

  names(x = transcripts_per_gene) <- c("Gene", "Transcripts")
  
  transcripts_per_gene$Tool <- tools[i]
  
  all_data_TX <- dplyr::bind_rows(all_data_TX, exons_per_transcript)
  all_data_GN <- dplyr::bind_rows(all_data_GN, transcripts_per_gene)
}

plot <- ggplot(data = all_data_TX, aes(x = Exons, fill = Tool)) +
  coord_cartesian(xlim = c(1, 10)) +
  facet_wrap(~Tool, scales = "free") +
  geom_histogram(binwidth = 1, colour = "black") +
  labs(x = "Exons per transcript", y = "Number of transcripts") +
  scale_fill_nejm() +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = scales::comma) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        text = element_text(family = "Helvetica"))

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.2 Exons per transcript/Exons per transcript.pdf", plot = plot, width = 11, height = 8.5)

wide_data <- all_data_GN %>% pivot_wider(names_from = Tool, values_from = Transcripts) %>% drop_na()
wide_data$rowSums <- rowSums(x = wide_data[, -1])

mat <- cor(x = wide_data[, 2:5])

pheatmap(mat = mat, 
         color = brewer.pal(n = 6, name = "Purples"), 
         main = "Number of known transcripts per gene",
         fontsize = 20,
         filename = "~/Desktop/Nanopore/Benchmarking/3.2 Exons per transcript/Correlation heatmap.pdf", 
         width = 11.5, 
         height = 8)