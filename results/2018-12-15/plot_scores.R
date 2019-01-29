library(ggplot2)

scores <- read.table('table_scores.txt', header=FALSE, col.names=c('numSamples', 'Sample', 'pca1', 'pca2', 'longitude', 'latitude'))
scores$pca1_norm <- numeric(length=length(scores$pca1))
#attach(scores)

for (i in min(scores$numSamples):max(scores$numSamples)) {
   filter <- scores$numSamples == i
   scores$pca1_norm[filter] <- (scores$pca1[filter] - min(scores$pca1[filter])) / (max(scores$pca1[filter]) - min(scores$pca1[filter]))
}

png(filename = 'scores.png', width=1000, height=1000)
g <- ggplot()
for (i in levels(scores$Sample)) {
   filter <- scores$Sample == i
   g <- g + geom_line(data = scores[filter,], mapping = aes(x = numSamples, y = pca1_norm, color = longitude))
}

g

dev.off()
#detach(scores)
write.table(scores, file = "table_scores_norm.txt", append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE, qmethod = c("escape", "double"),
                 fileEncoding = "")

