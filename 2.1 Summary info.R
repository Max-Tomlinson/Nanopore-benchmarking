library(dplyr)
library(ggplot2)
library(ggsci)
library(jsonlite)
library(RColorBrewer)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

file <- readLines(con = "~/Desktop/Nanopore/Benchmarking/2.1 Summary info/Summary_NanoStats.txt")
start <- grep(pattern = "^Reads", x = file)

quality <- file[(start)]
quality <- strsplit(x = quality, split = "\t")

data <- lapply(X = quality, FUN = function(x) {
  QScore      <- gsub(pattern = ".*>(Q\\d+):.*", replacement = "\\1", x = x[1])
  Counts      <- sapply(X = strsplit(x = x[-1], split = " "), FUN = function(d) d[1])
  Percentages <- sapply(X = strsplit(x = x[-1], split = " "), FUN = function(d) gsub(pattern = "[()%]", replacement = "", x = d[2]))
  Megabases   <- sapply(X = strsplit(x = x[-1], split = " "), FUN = function(d) gsub(pattern = "Mb", replacement = "", x = d[3]))
  list(QScore = QScore, Counts = Counts, Percentages = Percentages, Megabases = Megabases)
})

data <- do.call(what = rbind, args = lapply(X = data, FUN = function(x) {
  data.frame(QScore = x$QScore, Counts = x$Counts, Percentages = as.numeric(x = x$Percentages), Megabases = as.numeric(x = x$Megabases))
}))

machine <- readLines(con = "~/Desktop/Nanopore/Benchmarking/2.1 Summary info/Summary_NanoStats.txt", n = 1)
machine <- unlist(x = strsplit(x = machine, split = "\\s+"))
machine <- machine[-1]
machine <- sapply(X = strsplit(x = machine, split = "_"), FUN = function(x) x[3])

data$Machine <- rep(x = machine, 5)
data$Machine <- factor(x = data$Machine, levels = unique(x = data$Machine))
data$QScore <- factor(x = data$QScore, levels = c("Q15", "Q12", "Q10", "Q7", "Q5"))

plot <- ggplot(data = data, aes(x = Machine, y = Percentages, fill = QScore)) +
  geom_bar(stat = "identity", position = "identity", width = 0.7, color = "black") +
  labs(x = NULL, y =  "Percentage of reads above quality cut-off") +
  scale_fill_nejm() +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 16),
        legend.title = element_blank())

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.1 Summary info/Barplot percentage read quality.pdf", plot = plot, width = 11, height = 8.5, units = "in", useDingbats = FALSE)

OverlayLogHistogram <- fromJSON(txt = "~/Desktop/Nanopore/Benchmarking/2.1 Summary info/Summary_NanoComp_OverlayLogHistogram.json")
OverlayLogHistogram <- OverlayLogHistogram$data

data <- data.frame(x = numeric(), y = numeric(), element = factor())

for(i in seq_along(along.with = OverlayLogHistogram$x)) {
  temp <- data.frame(
    x = OverlayLogHistogram$x[[i]],
    y = OverlayLogHistogram$y[[i]],
    element = factor(x = sapply(X = strsplit(x = OverlayLogHistogram$name[i], split = "_"), FUN = function(x) x[3]))
  )
  data <- rbind(data, temp)
}

plot <- ggplot(data = data, aes(x = x, y = y, colour = "white", fill = element)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 0) +
  labs(x = "Read length", y = "Number of reads") +
  scale_x_continuous(labels = function(x) 10^x, breaks = log10(c(1, 10, 100, 1000, 10000))) +
  theme(legend.position = "none")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.1 Summary info/Barplot read lengths.pdf", plot = plot, width = 11, height = 8.5, units = "in", useDingbats = FALSE)

log_length_violin <- fromJSON(txt = "~/Desktop/Nanopore/Benchmarking/2.1 Summary info/Summary_NanoComp_log_length_violin.json")
log_length_violin <- log_length_violin$data

data <- data.frame(x = factor(), y = numeric())

for(i in seq_along(along.with = log_length_violin$x)) {
  temp <- data.frame(
    x = factor(sapply(X = strsplit(x = log_length_violin$name[i], split = "_"), function(x) x[3])),
    y = log_length_violin$y[[i]]
  )
  data <- rbind(data, temp)
}
 
plot <- ggplot(data = data, aes(x = x, y = y, fill = x)) +
  geom_violin(scale = "width") +
  geom_boxplot(outlier.colour = NA, width = 0.25, fill = "white") +
  labs(x = NULL, y = "Read length") +
  scale_fill_manual(values = brewer.pal(n = 9, name = "Set1")) +
  scale_y_continuous(labels = function(x) 10^x, breaks = log10(c(1, 10, 100, 1000, 10000))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.1 Summary info/Violin plot read lengths.pdf", plot = plot, width = 11, height = 8.5, units = "in", useDingbats = FALSE)

quals_violin <- fromJSON(txt = "~/Desktop/Nanopore/Benchmarking/2.1 Summary info/Summary_NanoComp_quals_violin.json")
quals_violin <- quals_violin$data

data <- data.frame(x = numeric(), y = numeric(), element = factor())

for(i in seq_along(along.with = quals_violin$x)) {
  temp <- data.frame(
    x = factor(x = sapply(X = strsplit(x = quals_violin$name[i], split = "_"), function(x) x[3])),
    y = quals_violin$y[[i]]
  )
  data <- rbind(data, temp)
}

plot <- ggplot(data = data, aes(x = x, y = y, fill = x)) +
  geom_violin(scale = "width") +
  geom_boxplot(outlier.colour = NA, width = 0.25, fill = "white") +
  labs(x = NULL, y = "Qscore") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(filename = "~/Desktop/Nanopore/Benchmarking/2.1 Summary info/Violin plot quality scores.pdf", plot = plot, width = 11, height = 8.5, units = "in", useDingbats = FALSE)