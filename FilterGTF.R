library(data.table)

cts <- fread(input = "~/Desktop/Nanopore/References/bambu/counts_transcript.txt")
gtf <- fread(input = "~/Desktop/Nanopore/References/bambu/extended_annotations.gtf")
gtf$transcript_id <- gsub(pattern = '.*transcript_id "([^"]+)";.*', replacement = '\\1', x = gtf$V9)

expressed_ids <- cts$TXNAME[rowSums(x = cts[, -c(1,2), with = FALSE]) > 0]
expressed_gtf <- gtf[gtf$transcript_id %in% expressed_ids, ]
expressed_gtf <- expressed_gtf[V7 %in% c("+", "-"), ]

fwrite(x = expressed_gtf[, 1:9], file = "~/Desktop/Nanopore/References/annotations_bambu.gtf", quote = FALSE, sep = "\t", col.names = FALSE)

cts <- fread(input = "~/Desktop/Nanopore/References/FLAIR/flair.quantify.counts.tsv")
cts$ids <- gsub(pattern = "_.*", replacement = "", x = cts$ids)
gtf <- fread(input = "~/Desktop/Nanopore/References/FLAIR/firstpass.isoforms.gtf")
gtf$transcript_id <- gsub(pattern = '.*transcript_id "([^"]+)";.*', replacement = '\\1', x = gtf$V9)

expressed_ids <- cts$ids[rowSums(x = cts[, -1, with = FALSE]) > 0]
expressed_gtf <- gtf[gtf$transcript_id %in% expressed_ids, ]
expressed_gtf <- expressed_gtf[V5 >= V4 - 1]
expressed_gtf <- expressed_gtf[V7 %in% c("+", "-"), ]
expressed_gtf[, V9 := gsub("\\|", "_", V9)]

fwrite(x = expressed_gtf[, 1:9], file = "~/Desktop/Nanopore/References/annotations_flair.gtf", quote = FALSE, sep = "\t", col.names = FALSE)

cts <- fread(input = "~/Desktop/Nanopore/References/isoquant/OUT.transcript_model_grouped_counts.tsv")
gtf <- fread(input = "~/Desktop/Nanopore/References/isoquant/OUT.transcript_models.gtf", header = FALSE)
gtf$transcript_id <- gsub(pattern = '.*transcript_id "([^"]+)";.*', replacement = '\\1', x = gtf$V9)
gtf <- gtf[V3 != "gene", ]

expressed_ids <- cts$`#feature_id`[rowSums(x = cts[, -1, with = FALSE]) > 0]
expressed_gtf <- gtf[gtf$transcript_id %in% expressed_ids,]
expressed_gtf <- expressed_gtf[V7 %in% c("+", "-"), ]

fwrite(x = expressed_gtf[, 1:9], file = "~/Desktop/Nanopore/References/annotations_isoquant.gtf", quote = FALSE, sep = "\t", col.names = FALSE)

cts <- fread(input = "~/Desktop/Nanopore/References/stringtie/transcript_count_matrix.csv")
gtf <- fread(input = "~/Desktop/Nanopore/References/stringtie/stringtie.merged.gtf", header = FALSE)
gtf$transcript_id <- gsub(pattern = '.*transcript_id "([^"]+)";.*', replacement = '\\1', x = gtf$V9)

expressed_ids <- cts$transcript_id[rowSums(x = cts[, -1, with = FALSE]) > 0]
expressed_gtf <- gtf[gtf$transcript_id %in% expressed_ids,]
expressed_gtf <- expressed_gtf[V7 %in% c("+", "-"), ]

fwrite(x = expressed_gtf[, 1:9], file = "~/Desktop/Nanopore/References/annotations_stringtie.gtf", quote = FALSE, sep = "\t", col.names = FALSE)