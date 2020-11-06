#
# Usage:
#
#   R --no-save --slave < fabbababa.R --args <input file | input dir> [ > output_file]
#
# If the only expected argument is a directory, all files named "rep*.freq"
# will be processed. If it is a file, only that one will be processed.
#
# ================================================================================
# FUNCTIONS
# ================================================================================

abba <- function(p1, p2, p3) (1 - p1) * p2 * p3

baba <- function(p1, p2, p3) p1 * (1 - p2) * p3

D.stat <- function(p1, p2, p3) {
   (sum(abba(p1, p2, p3)) - sum(baba(p1, p2, p3))) /
   (sum(abba(p1, p2, p3)) + sum(baba(p1, p2, p3)))
}

f_hom <- function(p1, p2, p3) {
   (sum(abba(p1, p2, p3)) - sum(baba(p1, p2, p3))) /
   (sum(abba(p1, p3, p3)) - sum(baba(p1, p3, p3)))
}

f_d <- function(p1, p2, p3) {
   (sum(abba(p1, p2, p3)) - sum(baba(p1, p2, p3))) /
   (sum(abba(p1, pmax(p2,p3), pmax(p2,p3))) - sum(baba(p1, pmax(p2,p3), pmax(p2,p3))))
}

# f_dM is defined in the supplementary materials of Malinsky et al. 2015, Science 350, 1493.
# Unfortunately, not very well specified, since it defines two f_dM, for different
# subsets of sites. I suppose the idea is to add them up.
f_dM <- function(p1, p2, p3) {
   z <- p2 > p1
   (sum(abba(p1,p2,p3)) - sum(baba(p1,p2,p3))) /
   ( sum(abba(p1[z], p2[z], p3[z])) - sum(baba(p1[z], pmax(p2[z],p3[z]), pmax(p2[z],p3[z]))) -
     sum(abba(pmax(p1[!z],p3[!z]), p2[!z], pmax(p1[!z],p3[!z]))) +
     sum(baba(pmax(p1[!z],p3[!z]), p2[!z], pmax(p1[!z],p3[!z]))) )
}

#==================================================================================
#
#==================================================================================

arguments <- commandArgs(trailingOnly = TRUE)
OutputTable <- array(numeric(), c(0,4))  # First dimension of size 0, to initialize it empty.
if (file.info(arguments[1])$isdir) {
   for (File in dir(arguments[1], pattern="rep.+freq")) {
#      if (length(system2('lsof', args=paste(arguments[1], File, sep='/'), stdout=TRUE)) == 0) {
         FT <- read.table(paste(arguments[1], File, sep="/"), header=TRUE,
            colClasses = c('NULL', 'NULL', 'numeric', 'numeric', 'numeric', 'NULL'))
         names(FT) <- c('P1', 'P2', 'P3')
         OutputTable <- rbind(OutputTable,
            c( D.stat(FT$P1, FT$P2, FT$P3),
               f_hom(FT$P1, FT$P2, FT$P3),
               f_d(FT$P1, FT$P2, FT$P3),
               f_dM(FT$P1, FT$P2, FT$P3) ))
#      }
   }
} else {
   FT <- read.table(arguments[1], header=TRUE,
      colClasses = c('NULL', 'NULL', 'numeric', 'numeric', 'numeric', 'NULL'))
   names(FT) <- c('P1', 'P2', 'P3')
   OutputTable <- rbind(OutputTable,
      c( D.stat(FT$P1, FT$P2, FT$P3),
         f_hom(FT$P1, FT$P2, FT$P3),
         f_d(FT$P1, FT$P2, FT$P3),
         f_dM(FT$P1, FT$P2, FT$P3) ))
}

write.table(as.matrix(OutputTable), col.names=FALSE, row.names=FALSE)

