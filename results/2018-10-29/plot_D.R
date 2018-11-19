library(ggplot2)

genes <- read.table('genes.abbababa.csv', header=TRUE, sep=',')
inter <- read.table('inter.abbababa.csv', header=TRUE, sep=',')

D <- data.frame(region=c(rep('genic', length(genes$D)), rep('intergenic', length(inter$D))),
                D=c(genes$D, inter$D))

png(filename='D.png')
ggplot(data=D, mapping=aes(x=region, y=D)) + geom_boxplot()
dev.off()
