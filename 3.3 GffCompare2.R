library(GenomeInfoDb)
library(ggplot2)
library(Gviz)
library(magrittr)
library(RColorBrewer)
library(rtracklayer)
library(scales)
library(tidyverse)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

tracking_files <- list.files(path = "~/Desktop/Nanopore/Benchmarking/3.1+3.3 GffCompare", pattern = "\\.tracking$", full.names = TRUE)[-3]
tracking <- data.frame()

for (file in tracking_files) {
  tool_name <- sub(pattern = ".*/(.*)\\.tracking$", replacement = "\\1", x = file)
  tool_name <- switch(EXPR = tool_name, "bambu" = "Bambu", "flair" = "FLAIR", "isoquant" = "IsoQuant", "stringtie" = "StringTie")
  
  tracking_file <- read_tsv(file = file, col_names = c("Query_transfrag_id", "Query_locus_id", "Reference_gene_id", "Class_code", "qJ"), show_col_types = FALSE)
  
  if (tool_name == "FLAIR") {
    tracking_file$qJ <- gsub(pattern = "\\|(\\w{8})\\b", replacement = "_\\1", x = tracking_file$qJ)
  }
  
  tracking_file <- separate(data = tracking_file, col = qJ, into = c("gene_id", "transcript_id", "num_exons", "FPKM", "TPM", "cov", "len"), sep = "\\|", convert = TRUE)
  tracking <- bind_rows(tracking, tracking_file %>% mutate(Tool = tool_name))
}

tracking <- tracking %>%
  group_by(Tool, Class_code) %>%
  filter(n() >= 2) %>%
  ungroup()

summary <- tracking %>% 
  group_by(Tool, Class_code) %>% 
  summarise(
    length = mean(len, na.rm = TRUE),
    log2 = mean(log2(x = len), na.rm = TRUE),
    sd_length = sd(len, na.rm = TRUE),
    exons = mean(num_exons, na.rm = TRUE),
    sd_exons = sd(num_exons, na.rm = TRUE),
  )

p1 <- ggplot(data = tracking, aes(x = log2(len), fill = Class_code)) +
  facet_wrap(~Tool) +
  geom_density(alpha = 0.5) +
  labs(x = "Log transcript length", y = "Density", fill = "Class code") +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme(legend.text = element_text(size = 16))

p2 <- ggplot(data = tracking, aes(x = num_exons, fill = Class_code)) +
  coord_cartesian(xlim = c(1, 10)) +
  facet_wrap(~Tool, scales = "free_y") +
  geom_density(alpha = 0.5) +
  labs(x = "Exons per transcript", y = "Density", fill = "Class code") +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme(legend.text = element_text(size = 16))

tracking$Annotation <- ifelse(test = startsWith(x = tracking$transcript_id, prefix = "ENST"), yes = "Known", no = "Novel")
tracking <- tracking %>%
  group_by(Tool, gene_id, Annotation) %>%
  summarise(Number = n(), .groups = "drop") %>%
  ungroup() %>%
  pivot_wider(names_from = "Annotation", values_from = "Number", values_fill = 0)
tracking$Novel <- ifelse(test = tracking$Novel == 0, yes = '0', no = ifelse(test = tracking$Novel == 1, yes = '1', no = ifelse(test = tracking$Novel == 2, yes = '2', no = "3+")))
tracking$Known <- ifelse(test = tracking$Known == 0, yes = '0', no = ifelse(test = tracking$Known == 1, yes = '1', no = ifelse(test = tracking$Known == 2, yes = '2', no = "3+")))

tracking <- tracking %>% 
  group_by(Tool, Novel, Known) %>% 
  summarise(sum = n()) %>% 
  ungroup() %>% 
  mutate(Novel = factor(x = Novel, levels = c("0", "1", "2", "3+"))) %>% 
  mutate(Known = factor(x = Known, levels = c("0", "1", "2", "3+")))

p3 <- ggplot(data = tracking, aes(x = Known, y = Novel, fill = log2(x = sum))) +
  facet_wrap(~Tool) +
  geom_tile() +
  labs(x = "Annotated transcripts per gene", y = "Novel transcripts per gene", fill = "Log") +
  scale_fill_gradient2() +
  theme(legend.position = "bottom")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.1+3.3 GffCompare/Combined lengths.pdf", plot = p1, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.1+3.3 GffCompare/Combined exons.pdf", plot = p2, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Benchmarking/3.1+3.3 GffCompare/Combined heatmap.pdf", plot = p3, width = 11, height = 8.5)

tool_names <- c("bambu", "flair", "isoquant", "stringtie")

pdf("~/Desktop/Nanopore/Benchmarking/3.1+3.3 GffCompare/Gene plots.pdf", width = 11, height = 8.5)
par(mfrow = c(1, length(x = tool_names)), mar = c(5, 4, 4, 2) + 0.1)

for (tool_name in tool_names) {
  
  gtf <- import.gff(con = paste0("~/Desktop/Nanopore/Benchmarking/3.1+3.3 GffCompare/", tool_name, ".annotated.gtf"))
  gtf <- as.data.frame(x = gtf)
 
  if (tool_name == "stringtie") {
    gtf$gene_id <- sub("\\..*$", "", gtf$ref_gene_id)
  } else {
    gtf$gene_id <- sub("\\..*$", "", gtf$gene_id)
  }
  
  gtf$label <- ifelse(test = grepl("^ENST", gtf$transcript_id), yes = "Known", no = "Novel")
  gtf <- gtf[which(x = gtf$gene_id == "ENSG00000160789"), ]
  
  gene_track <- GeneRegionTrack(range = gtf, name = tool_name, showId = FALSE, genome = "hg38", chromosome = "chr1")
  
  plotTracks(gene_track, from = min(gtf$start), to = max(gtf$end), main = "LMNA")
}

dev.off()
