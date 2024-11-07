library(Biostrings)
library(dplyr)
library(IRanges)
library(Repitools)
library(rtracklayer)
library(UpSetR)

gen <- readDNAStringSet("~/Desktop/Nanopore/References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
gen_len <- width(x = gen)
names(x = gen_len) <- sub(pattern = " .*", replacement = "", x = names(x = gen))
gen_bin <- genomeBlocks(genome = gen_len, width = 2000)

annotation <- function(file) {
  annot <- import(con = file) %>% .[.$type == "transcript", ]
  annot <- subsetByOverlaps(x = annot, ranges = gen_bin)
  annot <- annot[annot$type == "transcript", ]
  return(annot)
}

Bambu <- annotation(file = "~/Desktop/Nanopore/References/annotations_bambu.gtf")
FLAIR <- annotation(file = "~/Desktop/Nanopore/References/annotations_flair.gtf")
IsoQuant <- annotation(file = "~/Desktop/Nanopore/References/annotations_isoquant.gtf")
StringTie <- annotation(file = "~/Desktop/Nanopore/References/annotations_StringTie.gtf")

overlaps <- function(annotation) {
  overlap <- queryHits(x = findOverlaps(query = gen_bin, subject = annotation))
  return(overlap)
}

Bambu <- overlaps(annotation = Bambu)
FLAIR <- overlaps(annotation = FLAIR)
IsoQuant <- overlaps(annotation = IsoQuant)
StringTie <- overlaps(annotation = StringTie)

overlap <- list(Bambu = Bambu, FLAIR = FLAIR, IsoQuant = IsoQuant, StringTie = StringTie)

bin_mat <- fromList(input = overlap)

pdf(file = "~/Desktop/Nanopore/Benchmarking/3.4 Upset.pdf", width = 11, height = 8.5,  family = "Helvetica")
options(scipen = 999)

upset(
  data = bin_mat,
  mainbar.y.label = "Overlapping isoform coordinates",
  point.size = 4,
  order.by = "freq",
  text.scale = c(3, 2.5, 2.5, 1.5, 2.5, 1.5)
)

dev.off()