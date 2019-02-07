# I use the accurate statistics obtained from non-overlapping windows. Unfortunately,
# there are no abba-baba parameter estimates in the same windows.
popgen     <- read.table('erin63.AccuStats.csv',  header=TRUE, sep=',')
geneStats  <- read.table('genes.PopGenStats.csv', header=TRUE, sep=',')
#geneABBA   <- read.table('genes.abbababa.csv',    header=TRUE, sep=',')  # not used
interStats <- read.table('inter.PopGenStats.csv', header=TRUE, sep=',')
#interABBA  <- read.table('inter.abbababa.csv',    header=TRUE, sep=',')  # not used

library(ggplot2)
library(gridExtra)

p1 <- ggplot(data=popgen, mapping=aes(x=pi_roumanicus, y=dxy_roumanicus_europaeus)) + 
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. roumanicus') +
      ylab('Divergence')

p2 <- ggplot(data=geneStats, mapping=aes(x=pi_roumanicus, y=dxy_roumanicus_europaeus)) + 
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. roumanicus') +
      ylab('Divergence') + ggtitle('Genic regions')

p3 <- ggplot(data=interStats, mapping=aes(x=pi_roumanicus, y=dxy_roumanicus_europaeus)) +
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. roumanicus') +
      ylab('Divergence') + ggtitle('Intergenic regions')
  
p4 <- ggplot(data=popgen, mapping=aes(x=pi_europaeus, y=dxy_roumanicus_europaeus)) + 
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. europaeus') +
      ylab('Divergence')     

p5 <- ggplot(data=geneStats, mapping=aes(x=pi_europaeus, y=dxy_roumanicus_europaeus)) +
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. europaeus') +
      ylab('Divergence') + ggtitle('Genic regions')

p6 <- ggplot(data=interStats, mapping=aes(x=pi_europaeus, y=dxy_roumanicus_europaeus)) +
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. europaeus') +
      ylab('Divergence') + ggtitle('Intergenic regions')

png(filename='pi_dxy.png')
grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2)
dev.off()

