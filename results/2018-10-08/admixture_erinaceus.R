k2 <- as.matrix(read.table('erinaceus.2.Q'))
k3 <- as.matrix(read.table('erinaceus.3.Q'))
k4 <- as.matrix(read.table('erinaceus.4.Q'))
k5 <- as.matrix(read.table('erinaceus.5.Q'))
k6 <- as.matrix(read.table('erinaceus.6.Q'))
k7 <- as.matrix(read.table('erinaceus.7.Q'))
k8 <- as.matrix(read.table('erinaceus.8.Q'))
k9 <- as.matrix(read.table('erinaceus.9.Q'))

names <- as.matrix(read.table('out.012.indv'))[,1]

library('RColorBrewer')

# This gives the order of individuals according to decreasing proportions
# of ancestry 1 (europaeus), 2 (roumanicus), 3 (concolor) from k=3, and then
# according to decreasing proportions of ancestry 4 (Western europaeus), 1
# (Eastern europaeus), 3 (roumanicus), and 2 (concolor) from k=4.

z1 <- cbind(k3, k4)
ordered <- order(z1[,1], z1[,2], z1[,3], z1[,7], z1[,4], z1[,6], z1[,5], decreasing=TRUE)

# To get the labels vertical, we need to pass the las=2 option to par. Otherwise,
# barplot's automatic settings override the option.


barX  <- barplot(t(k2[ordered,]), col=brewer.pal(12,'Set3'), ylab="Ancestry", border=NA, cex.names=0.8, main="K=2")
text((barX-0.7),-0.12,cex=0.7, names[ordered], xpd=TRUE, srt=45)

barX <- barplot(t(k3[ordered,]), col=brewer.pal(12,'Set3'), ylab="Ancestry", border=NA, cex.names=0.8, main="K=3")
text((barX-0.7),-0.12,cex=0.7, names[ordered], xpd=TRUE, srt=45)

barX <- barplot(t(k4[ordered,]), col=brewer.pal(12,'Set3'), ylab="Ancestry", border=NA, cex.names=0.8, main="K=4")
text((barX-0.7),-0.12,cex=0.7, names[ordered], xpd=TRUE, srt=45)

barX <- barplot(t(k5[ordered,]), col=brewer.pal(12,'Set3'), ylab="Ancestry", border=NA, cex.names=0.8, main="K=5")
text((barX-0.7),-0.12,cex=0.7, names[ordered], xpd=TRUE, srt=45)

barX <- barplot(t(k6[ordered,]), col=brewer.pal(12,'Set3'), ylab="Ancestry", border=NA, cex.names=0.8, main="K=6")
text((barX-0.7),-0.12,cex=0.7, names[ordered], xpd=TRUE, srt=45)

barX <- barplot(t(k7[ordered,]), col=brewer.pal(12,'Set3'), ylab="Ancestry", border=NA, cex.names=0.8, main="K=7")
text((barX-0.7),-0.12,cex=0.7, names[ordered], xpd=TRUE, srt=45)

barX <- barplot(t(k8[ordered,]), col=brewer.pal(12,'Set3'), ylab="Ancestry", border=NA, cex.names=0.8, main="K=8")
text((barX-0.7),-0.12,cex=0.7, names[ordered], xpd=TRUE, srt=45)

barX <- barplot(t(k9[ordered,]), col=brewer.pal(12,'Set3'), ylab="Ancestry", border=NA, cex.names=0.7, main="K=9")
text((barX-0.7),-0.12,cex=0.7, names[ordered], xpd=TRUE, srt=45)

dev.off()
