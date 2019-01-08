# Usage:
#
# R [-q] [--no-save|--save] <run_Adegenet.R --args <012_file> <num_sites> <indv_file> <positions_file> <output_file>

library('ape')
library('pegas')
library('seqinr')
library('ggplot2')
library('adegenet')

arguments <- commandArgs(trailingOnly=TRUE)

SNP <- read.table(file=arguments[1], header=FALSE, sep="\t",
   colClasses=c('NULL', rep('integer', arguments[2])))
SNPgl <- new('genlight', SNP)
indNames(SNPgl) <- as.matrix(read.table(arguments[3]))
Positions <- read.table(arguments[4], header=FALSE,
   col.names=c('chr', 'pos'))
chromosome(SNPgl) <- Positions$chr
position(SNPgl) <- Positions$pos
pca1 <- glPca(SNPgl, nf=2)
write.table(pca1$scores, file=arguments[5], row.names=TRUE, quote=FALSE, col.names=FALSE)
