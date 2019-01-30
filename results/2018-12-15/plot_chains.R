library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly=TRUE)
infile1 <- args[1]
infile2 <- args[2]
outfile <- args[3]

chain1 <- data.frame(t(as.matrix(read.table(infile1, sep=','))))
chain2 <- data.frame(t(as.matrix(read.table(infile2, sep=','))))
chain1$pos <- 1:dim(chain1)[1]
chain2$pos <- 1:dim(chain2)[1]

# If the parameter is LnL, the input file has one row, which should become one column in
# the data frame. But if its alpha or beta, it has as many rows as loci.

if (dim(chain1)[2] >= 9) {
   subsetCols <- sample(dim(chain1)[2], size = 9)
   p1 <- ggplot() + geom_point(data=chain1, mapping=aes(x=chain1$pos, y=chain1[,subsetCols[1]]), colour='red', alpha=0.5) + 
                    geom_point(data=chain2, mapping=aes(x=chain2$pos, y=chain2[,subsetCols[1]]), colour='blue', alpha=0.4)
   p2 <- ggplot() + geom_point(data=chain1, mapping=aes(x=chain1$pos, y=chain1[,subsetCols[2]]), colour='red', alpha=0.5) + 
                    geom_point(data=chain2, mapping=aes(x=chain2$pos, y=chain2[,subsetCols[2]]), colour='blue', alpha=0.4)
   p3 <- ggplot() + geom_point(data=chain1, mapping=aes(x=chain1$pos, y=chain1[,subsetCols[3]]), colour='red', alpha=0.5) + 
                    geom_point(data=chain2, mapping=aes(x=chain2$pos, y=chain2[,subsetCols[3]]), colour='blue', alpha=0.4)
   p4 <- ggplot() + geom_point(data=chain1, mapping=aes(x=chain1$pos, y=chain1[,subsetCols[4]]), colour='red', alpha=0.5) + 
                    geom_point(data=chain2, mapping=aes(x=chain2$pos, y=chain2[,subsetCols[4]]), colour='blue', alpha=0.4)
   p5 <- ggplot() + geom_point(data=chain1, mapping=aes(x=chain1$pos, y=chain1[,subsetCols[5]]), colour='red', alpha=0.5) + 
                    geom_point(data=chain2, mapping=aes(x=chain2$pos, y=chain2[,subsetCols[5]]), colour='blue', alpha=0.4)
   p6 <- ggplot() + geom_point(data=chain1, mapping=aes(x=chain1$pos, y=chain1[,subsetCols[6]]), colour='red', alpha=0.5) + 
                    geom_point(data=chain2, mapping=aes(x=chain2$pos, y=chain2[,subsetCols[6]]), colour='blue', alpha=0.4)
   p7 <- ggplot() + geom_point(data=chain1, mapping=aes(x=chain1$pos, y=chain1[,subsetCols[7]]), colour='red', alpha=0.5) + 
                    geom_point(data=chain2, mapping=aes(x=chain2$pos, y=chain2[,subsetCols[7]]), colour='blue', alpha=0.4)
   p8 <- ggplot() + geom_point(data=chain1, mapping=aes(x=chain1$pos, y=chain1[,subsetCols[8]]), colour='red', alpha=0.5) + 
                    geom_point(data=chain2, mapping=aes(x=chain2$pos, y=chain2[,subsetCols[8]]), colour='blue', alpha=0.4)
   p9 <- ggplot() + geom_point(data=chain1, mapping=aes(x=chain1$pos, y=chain1[,subsetCols[9]]), colour='red', alpha=0.5) + 
                    geom_point(data=chain2, mapping=aes(x=chain2$pos, y=chain2[,subsetCols[9]]), colour='blue', alpha=0.4)
   g <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow=3)
} else {
   if (dim(chain1)[2] == 2) {
      g <- ggplot() + geom_point(data=chain1, mapping=aes(x=chain1$pos, y=chain1[,1]), colour='red', alpha=0.5) +
                      geom_point(data=chain2, mapping=aes(x=chain2$pos, y=chain2[,1]), colour='blue', alpha=0.4)
   }
}

ggsave(outfile, plot = g)
