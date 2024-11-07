library(cowplot)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggsci)
library(patchwork)
library(scales)
library(tidyverse)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

tools <- c("Bambu", "FLAIR", "IsoQuant", "StringTie")

Class_list <- list()
Count_list <- list()

for (i in 1:length(tools)) {
  Class_list[[tools[i]]] <- read.delim(file = list.files(path = "~/Desktop/Nanopore/Benchmarking/4.1 SQANTI3", pattern = "classification.txt", recursive = TRUE, full.names = TRUE)[i])
}

Class_list[["Bambu"]]$Method <- "Bambu"
Class_list[["FLAIR"]]$Method <- "FLAIR"
Class_list[["IsoQuant"]]$Method <- "IsoQuant"
Class_list[["StringTie"]]$Method <- "StringTie"

Count_list[["Bambu"]] <- read.delim(file = "~/Desktop/Nanopore/References/bambu/counts_transcript.txt")
Count_list[["Bambu"]]$Total <- rowSums(x = Count_list[["Bambu"]][, -c(1:2)], na.rm = TRUE)

Count_list[["FLAIR"]] <- read.delim(file = "~/Desktop/Nanopore/References/FLAIR/flair.quantify.counts.tsv")
Count_list[["FLAIR"]]$Total <- rowSums(x = Count_list[["FLAIR"]][, -1], na.rm = TRUE)
Count_list[["FLAIR"]]$ids <- gsub(pattern = "_.*", replacement = "", x = Count_list[["FLAIR"]]$ids)
Count_list[["FLAIR"]]$ids <- gsub(pattern = "\\|", replacement = "_", x = Count_list[["FLAIR"]]$ids)

Count_list[["IsoQuant"]] <- read.delim(file = "~/Desktop/Nanopore/References/IsoQuant/OUT.transcript_model_grouped_counts.tsv")
Count_list[["IsoQuant"]]$Total <- rowSums(x = Count_list[["IsoQuant"]][, -1], na.rm = TRUE)

Count_list[["StringTie"]] <- read.delim(file = "~/Desktop/Nanopore/References/stringtie/counts_transcript.txt", skip = 1)
Count_list[["StringTie"]]$Total <- rowSums(x = Count_list[["StringTie"]][, -c(1:8)], na.rm = TRUE)

Class_list[["Bambu"]]$Count <- Count_list[["Bambu"]][match(
  Class_list[["Bambu"]]$isoform,
  Count_list[["Bambu"]]$TXNAME
), "Total"]

Class_list[["FLAIR"]]$Count <- Count_list[["FLAIR"]][match(
  Class_list[["FLAIR"]]$isoform,
  Count_list[["FLAIR"]]$ids
), "Total"]

Class_list[["IsoQuant"]]$Count <- Count_list[["IsoQuant"]][match(
  Class_list[["IsoQuant"]]$isoform,
  Count_list[["IsoQuant"]]$X.feature_id
), "Total"]

Class_list[["StringTie"]]$Count <- Count_list[["StringTie"]][match(
  Class_list[["StringTie"]]$isoform,
  Count_list[["StringTie"]]$Geneid
), "Total"]

Class <- bind_rows(Class_list) %>%
  mutate(structural_category = recode(
    structural_category,
    "genic"        = "genic/intergenic",
    "genic_intron" = "genic/intergenic",
    "intergenic"   = "genic/intergenic"
  ))

data <- Class %>%
  mutate(structural_category = gsub(pattern = "_", replacement = " ", x = structural_category), structural_category = gsub(pattern = "^(.)", replacement = "\\U\\1", x = structural_category, perl = TRUE))

plot_class <- ggplot(data = data, aes(x = Method, fill = structural_category)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  labs(x = NULL, y = "Number of transcripts", fill = "Structural category") +
  scale_fill_nejm() +
  scale_y_continuous(labels = scales::comma)

data <- data %>% dplyr::count(Method, structural_category, wt = Count, name = "Count")

plot_count <- ggplot(data = data, aes(x = Method, y = Count, fill = structural_category)) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
  labs(x = NULL, y = "Number of reads", fill = "Structural category") +
  scale_fill_nejm() +
  scale_x_discrete() +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))

plot <- (plot_class + plot_count) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom", 
        legend.title = element_blank()) 

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/4.1 SQANTI3/Combined.pdf", plot = plot, width = 22, height = 8.5)

results_list <- list()

for (i in tools) {
  class_data <- Class[Class$Method == i, ]
  exonvector <- numeric(length = 5)
  intergenic <- nrow(x = filter(.data = class_data, structural_category == "intergenic"))
  
  for (n_exons in 1:4) {
    known <- nrow(x = filter(.data = class_data, grepl(pattern = "^ENST", x = isoform), exons == n_exons))
    novel <- nrow(x = filter(.data = class_data, !grepl(pattern = "^ENST", x = isoform), exons == n_exons))
    intergenic_known <- nrow(x = filter(.data = class_data, grepl(pattern = "^ENST", x = isoform), exons == n_exons, structural_category == "intergenic"))
    intergenic_novel <- nrow(x = filter(.data = class_data, !grepl(pattern = "^ENST", x = isoform), exons == n_exons, structural_category == "intergenic"))
    exonvector[n_exons] <- intergenic_novel
  }
  
  for (n_exons in 5) {
    intergenic_novel <- nrow(x = filter(.data = class_data, !grepl(pattern = "^ENST", x = isoform), exons >= n_exons, structural_category == "intergenic"))
    exonvector[n_exons] <- intergenic_novel
  }
  
  results_list[[i]] <- c(exonvector, intergenic)
}

results <- do.call(what = rbind, args = results_list)
rownames(x = results) <- names(x = Class_list)
colnames(x = results) <- c(paste0("Exon", 1:4), "Exon5+", "Total")

write.csv(x = results, "~/Desktop/Nanopore/Benchmarking/4.1 SQANTI3/Intergenic.csv", row.names = TRUE)

flair <- Class %>%
  filter(structural_category == 'genic/intergenic', Method == 'FLAIR', exons == 1, CDS_genomic_start > 0) #%>% 
  #select(where(~ all(complete.cases(.))))

flair_ranges <- GRanges(
  seqnames = Rle(values = flair$chrom),
  ranges = IRanges(start = pmin(flair$CDS_genomic_start, flair$CDS_genomic_end), end = pmax(flair$CDS_genomic_start, flair$CDS_genomic_end)),
  strand = Rle(values = flair$strand)
)

bambu <- Class %>%
  filter(Method == 'Bambu', CDS_genomic_start > 0) #%>% 
  #select(where(~ all(complete.cases(.))))

bambu_ranges <- GRanges(
  seqnames = Rle(values = bambu$chrom),
  ranges = IRanges(start = pmin(bambu$CDS_genomic_start, bambu$CDS_genomic_end), end = pmax(bambu$CDS_genomic_start, bambu$CDS_genomic_end)),
  strand = Rle(values = bambu$strand)
)

stringtie <- Class %>%
  filter(Method == 'StringTie', CDS_genomic_start > 0) #%>% 
  #select(where(~ all(complete.cases(.))))

stringtie_ranges <- GRanges(
  seqnames = Rle(values = stringtie$chrom),
  ranges = IRanges(start = pmin(stringtie$CDS_genomic_start, stringtie$CDS_genomic_end), end = pmax(stringtie$CDS_genomic_start, stringtie$CDS_genomic_end)),
  strand = Rle(values = stringtie$strand)
)

overlap1 <- findOverlaps(query = flair_ranges, subject = bambu_ranges, type = "equal")
overlap1 <- queryHits(x = overlap1)
overlap1 <- bambu[overlap1, ] %>% distinct()
overlap1$Overlap <- "Equal"

overlap2 <- findOverlaps(query = flair_ranges, subject = bambu_ranges, type = "within")
overlap2 <- queryHits(x = overlap2)
overlap2 <- bambu[overlap2, ] %>% distinct()
overlap2$Overlap <- "Within"

overlap3 <- findOverlaps(query = flair_ranges, subject = stringtie_ranges, type = "equal")
overlap3 <- queryHits(x = overlap3)
overlap3 <- stringtie[overlap3, ] %>% distinct()
overlap3$Overlap <- "Equal"

overlap4 <- findOverlaps(query = flair_ranges, subject = stringtie_ranges, type = "within")
overlap4 <- queryHits(x = overlap4)
overlap4 <- stringtie[overlap4, ] %>% distinct()
overlap4$Overlap <- "Within"

overlaps <- rbind(overlap1, overlap2, overlap3, overlap4)

data <- overlaps %>%
  group_by(Method, Overlap) %>%
  summarise(
    Transcript = sum(grepl(pattern = "^ENST", x = associated_transcript)) / n() * 100,
    Gene = sum(grepl(pattern = "^ENSG", x = associated_gene)) / n() * 100,
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Transcript, Gene), names_to = "Category", values_to = "Percentage")

plot <- ggplot(data = data, aes(x = Overlap, y = Percentage, fill = Category)) +
  facet_wrap(facets = ~Method) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") +
  labs(x = NULL, y = "Percentage of annotated features", fill = NULL) +
  scale_fill_nejm() +
  theme(legend.position = "bottom")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/4.1 SQANTI3/Overlap_annotated.pdf", plot = plot, width = 11, height = 8.5)

plot <- ggplot(data = overlaps[overlaps$Overlap == "Within", ], aes(x = exons, fill = Method)) +
  coord_cartesian(xlim = c(1, 10)) +
  geom_density(adjust = 2, alpha = 0.7, colour = "black") +
  labs(x = "Exons per transcript", y = "Density") +
  scale_fill_nejm() +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "none")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/4.1 SQANTI3/Overlap_exons.pdf", plot = plot, width = 11, height = 8.5)

data <- overlaps %>%
  mutate(structural_category = gsub(pattern = "_", replacement = " ", x = structural_category),
         structural_category = gsub(pattern = "^(.)", replacement = "\\U\\1", x = structural_category, perl = TRUE)) %>%
  mutate(subcategory = gsub(pattern = "_", replacement = " ", x = subcategory),
         subcategory = gsub(pattern = "^(.)", replacement = "\\U\\1", x = subcategory, perl = TRUE))

plot_class <- ggplot(data = data[data$Overlap == "Within", ], aes(x = Method, fill = structural_category)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  labs(x = NULL, y = "Number of transcripts", fill = "Structural category") +
  scale_fill_nejm() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "none")

data <- data %>%
  group_by(Method, Overlap, structural_category) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")

plot_count <- ggplot(data = data, aes(x = Overlap, y = Count, fill = structural_category)) +
  facet_wrap(facets = ~Method) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
  labs(x = NULL, y = "Number of reads", fill = "Structural category") +
  scale_fill_nejm() +
  scale_x_discrete() +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/4.1 SQANTI3/Overlap_combined.pdf", plot = (plot + plot_class), width = 11, height = 8.5)