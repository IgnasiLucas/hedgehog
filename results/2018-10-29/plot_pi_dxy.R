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
      ylab('Divergence') + ggtitle('A. Whole genome') + theme(plot.title=element_text(hjust=0.5))

p2 <- ggplot(data=geneStats, mapping=aes(x=pi_roumanicus, y=dxy_roumanicus_europaeus)) + 
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. roumanicus') +
      ylab('Divergence') + ggtitle('B. Genic regions') + theme(plot.title=element_text(hjust=0.5))

p3 <- ggplot(data=interStats, mapping=aes(x=pi_roumanicus, y=dxy_roumanicus_europaeus)) +
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. roumanicus') +
      ylab('Divergence') + ggtitle('C. Intergenic regions') + theme(plot.title=element_text(hjust=0.5))
  
p4 <- ggplot(data=popgen, mapping=aes(x=pi_europaeus, y=dxy_roumanicus_europaeus)) + 
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. europaeus') +
      ylab('Divergence')     

p5 <- ggplot(data=geneStats, mapping=aes(x=pi_europaeus, y=dxy_roumanicus_europaeus)) +
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. europaeus') +
      ylab('Divergence')

p6 <- ggplot(data=interStats, mapping=aes(x=pi_europaeus, y=dxy_roumanicus_europaeus)) +
      geom_point() + geom_smooth(method='lm') + xlab('Diversity in E. europaeus') +
      ylab('Divergence')

p <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2, widths=c(1.2, 1, 1))
ggsave('pi_dxy.pdf', plot=p, device='pdf', width=12, height=7, units='in')

# The figure saved above uses accurate statistics from non-overlapping windows, and
# distinguishing genic from intergenic regions. The two first plots would be more
# informative if points were color-coded to indicate the amount of introgression in
# the window where the statistics (diversity and divergence) was calculated. Unfortunately,
# the abba-baba statistics about introgression have not been calculated in the same
# windows, because of different requirements. I will repeat the plot above, now using the
# alternative data set where genic and intergenic regions are not distinguished, but
# introgression statistics are calculated in common windows with diversity and divergence,
# for the first, genome-wide plots. This has already been plotted separately in plot_pi_dxy_fd.R.

popgen     <- read.table('nonoverlap.PopGenStats.csv',  header=TRUE, sep=',',
              col.names=c('scaffold', 'start', 'end', 'mid', 'sites', 'pi_roumanicus',
              'pi_europaeus', 'dxy', 'Fst'))
abbababa   <- read.table('nonoverlap.abbababa.csv', header=TRUE, sep=',')
abbababa$fd[abbababa$D <= 0] <- 0.0
abbababa$fdM[abbababa$D <= 0] <- 0.0
allStats   <- merge(popgen, abbababa, by=c('scaffold', 'start', 'end'))
filter  <- !is.na(allStats$fd)

p1 <- ggplot(data=allStats[filter,], mapping=aes(x=pi_roumanicus, y=dxy)) +
      geom_point(mapping=aes(colour=fd), size=2) + geom_smooth(method='lm') +
      xlab('Diversity in E. roumanicus') + ylab('Divergence') +
      scale_colour_gradient2(low='yellow',mid='orange',high='red',midpoint=0.1) +
      ggtitle('A. Whole genome') + theme(plot.title=element_text(hjust=0.5))

p4 <- ggplot(data=allStats[filter,], mapping=aes(x=pi_europaeus, y=dxy)) +
      geom_point(mapping=aes(colour=fd), size=2) + geom_smooth(method='lm') +
      xlab('Diversity in E. europaeus') + ylab('Divergence') +
      scale_colour_gradient2(low='yellow',mid='orange',high='red',midpoint=0.1)

p <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2, widths=c(1.2, 1, 1))
ggsave('pi_dxy_2.pdf', plot=p, device='pdf', width=12, height=7, units='in')
