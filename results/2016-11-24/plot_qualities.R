z1 <- read.table('z1', col.names=c('qual','type'))
attach(z1)
png(filename='qualities.png', width=2000, height=500)
   boxplot(log(qual) ~ type, ylab="log(quality)", xlab="Type of variant")
dev.off()
detach(z1)

