library(caret)
library(DeconCell)
library(dplyr)
library(lubridate)
library(matrixStats)
library(tidyr)
library(truncnorm)

count.table <- read.delim(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Raw Counts/swap.merged.blood.gene.count.txt", row.names = 4, check.names = FALSE)[, -c(1:5)]

rownames(x = count.table) <- sub(pattern = "\\..*", replacement = "", x = rownames(x = count.table))

data("dCell.models")

dCell <- dCell.expProcessing(count.table = count.table)
dCell <- dCell.predict(dCell.exp = dCell, dCell.models = dCell.models)
dCell <- as.data.frame(x = dCell$dCell.prediction)

pca <- prcomp(x = dCell, center = TRUE, scale. = TRUE)
pca <- as.data.frame(x = pca$x)

colData <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/MetabPhenos_ALLEBSubjects_NoPhenoExclusions_WithHOMAs_20181120.txt", header = TRUE)
colData <- merge(x = colData, y = pca[, 1:3], by.x = 1, by.y = 0)

write.table(x = colData, file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/MetabPhenos_ALLEBSubjects_NoPhenoExclusions_WithHOMAs_20181120_dCell.txt", quote = FALSE, sep = "\t", row.names = FALSE)