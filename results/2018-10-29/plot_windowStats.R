library(ggplot2)
abbababa <- read.table('creh1.abbababa.csv',     header=TRUE, sep=',')
  popgen <- read.table('erin63.PopGenStats.csv', header=TRUE, sep=',')

# Windows where there is no excess of ABBA, the 'excess' cannot be quantified.
abbababa$fd[abbababa$D <= 0] <- 0.0

allStats <- merge(popgen, abbababa, by=c('scaffold', 'start', 'end'))

# I use only scaffolds with at least 27 data points.
goodScaff <- names(table(allStats$scaffold)[table(allStats$scaffold) >= 27])
filter <- allStats$scaffold %in% goodScaff

png(filename='windowStats.png', width=1000, height=1000)
ggplot(data=allStats[filter,]) + geom_line(mapping=aes(x=mid.x, y=dxy_roumanicus_europaeus), color='green') + 
   geom_line(mapping=aes(x=mid.x, y=pi_roumanicus), color='red') +
   geom_line(mapping=aes(x=mid.x, y=pi_europaeus),  color='blue') +
   geom_line(mapping=aes(x=mid.y, y=fd)) +
   facet_wrap(~scaffold) + xlab('position') + ylab('')
dev.off()
