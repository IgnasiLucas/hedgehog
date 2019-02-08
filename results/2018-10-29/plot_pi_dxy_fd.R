# Here I use statistics computed in common windows
popgen     <- read.table('nonoverlap.PopGenStats.csv',  header=TRUE, sep=',', 
              col.names=c('scaffold', 'start', 'end', 'mid', 'sites', 'pi_roumanicus', 'pi_europaeus', 'dxy', 'Fst'))
abbababa   <- read.table('nonoverlap.abbababa.csv', header=TRUE, sep=',')
abbababa$fd[abbababa$D <= 0] <- 0.0
abbababa$fdM[abbababa$D <= 0] <- 0.0
allStats   <- merge(popgen, abbababa, by=c('scaffold', 'start', 'end'))
filter  <- !is.na(allStats$fd)
filter2 <- !is.na(allStats$fd) & allStats$fd > 0
library(ggplot2)
library(gridExtra)

p1 <- ggplot(data=allStats[filter,], mapping=aes(x=pi_roumanicus, y=dxy)) + 
      geom_point(mapping=aes(colour=fd), size=2) + geom_smooth(method='lm') + xlab('Diversity in E. roumanicus') +
      ylab('Divergence') + scale_colour_gradient2(low='yellow',mid='orange',high='red',midpoint=0.1)

p2 <- ggplot(data=allStats[filter,], mapping=aes(x=pi_europaeus, y=dxy)) + 
      geom_point(mapping=aes(colour=fd), size=2) + geom_smooth(method='lm') + xlab('Diversity in E. europaeus') +
      ylab('Divergence') + scale_colour_gradient2(low='yellow',mid='orange',high='red',midpoint=0.1)

p3 <- ggplot(data=allStats[filter2,], mapping=aes(x=fd, y=pi_roumanicus)) + geom_point() + scale_x_continuous(trans='log10') + geom_smooth(method='lm')

p4 <- ggplot(data=allStats[filter2,], mapping=aes(x=fd, y=pi_europaeus))  + geom_point() + scale_x_continuous(trans='log10') + geom_smooth(method='lm')

p5 <- ggplot(data=allStats[filter2,], mapping=aes(x=fd, y=dxy)) + geom_point() + scale_x_continuous(trans='log10') + geom_smooth(method='lm')

p6 <- ggplot(data=allStats[filter,], mapping=aes(x=pi_roumanicus, y=pi_europaeus, colour=dxy)) + 
      geom_point(size=0.8) + scale_colour_gradient2(low='yellow',mid='orange',high='red',midpoint=0.45) +
      geom_smooth(method='lm')

png(filename='pi_dxy_fd.png')
grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3)
dev.off()

ggsave('pi_dxy_roumanicus.png', p1)
ggsave('pi_dxy_europaeus.png', p2)
