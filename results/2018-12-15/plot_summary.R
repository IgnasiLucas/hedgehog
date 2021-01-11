library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
print(args[3])
estimates <- read.table(args[1], sep=',', col.names=c('meanPar', 'medianPar', 'lowCI', 'highCI'))
positions <- read.table(args[2], col.names=c('contig', 'pos'))

if (dim(estimates)[1] == dim(positions)[1]) {
   estimates$contig <- positions$contig
   estimates$pos    <- positions$pos
   minSNPs <- 4
   goodScaff <- names(table(positions$contig)[table(positions$contig) >= minSNPs])
   while (length(goodScaff) > 25) {
      minSNPs <- minSNPs + 1
      goodScaff <- names(table(positions$contig)[table(positions$contig) >= minSNPs])
   }
   filter <- positions$contig %in% goodScaff

   g <- ggplot(data=estimates[filter,], mapping=aes(x=pos)) + geom_line(mapping=aes(y=meanPar), color='black') + 
           geom_line(mapping=aes(y=medianPar), color='red') +
           geom_ribbon(mapping=aes(ymin=lowCI, ymax=highCI), fill='gray', alpha=0.5) +
           facet_wrap(~contig) + xlab('position') + ylab('')
   ggsave(args[3], plot=g)
}
