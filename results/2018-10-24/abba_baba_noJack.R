# This version of the abba/baba test, taken originally from Simon Martin, does not run the jackknive
# process to estimate the standard error of the D statistic.
#
# =====================================================
#  USAGE
# =====================================================
#
# R [options] <abba_baba.R >outfile --args <freq_table> [<P1> <P2> <P3>]
#
#    options: --save, --no-save or --vanilla. One required.
#
#    freq_table  Tab-separated values of derived allele frequencies in populations P1, P2 and P3
#                It must have a header with the population names. Populations P1, P2, and P3
#                are asumed to correspond to columns 3, 4, and 6, after 'scaffold' and 'position'.
#
#    P1 P2 P3    If the population names are given as additional arguments, the order in which
#                they are specified will define which columns in freq_table correspond to P1, P2
#                and P3, as long as they coincide with the names in the header of freq_table.
#
#
# =====================================================
#  READING ARGUMENTS
# =====================================================
#
# The reason to use the syntax described above and commandArgs(), below, is to be able to
# use the same script for any set of populations, so that we don't need to edit the script
# every time we change the table of frequencies or the names of populations.

arguments <- commandArgs(trailingOnly=TRUE)

if (length(arguments) == 1) {
   freq_table <- read.table(arguments[1], header=TRUE, as.is=TRUE)
   if (dim(freq_table)[2] >= 5) {
      P1 <- names(freq_table)[3]
      P2 <- names(freq_table)[4]
      P3 <- names(freq_table)[5]
   } else {
      print("Bad format of freq_table. Scaffold and position required.")
      q(save="no", status=1)
   }
} else {
   if (length(arguments) == 4) {
      freq_table <- read.table(arguments[1], header=TRUE, as.is=TRUE)
      if (is.element(arguments[2], names(freq_table))) P1 <- arguments[2] else {paste(arguments[2], "not found in header."); q(save="no", status=1)}
      if (is.element(arguments[3], names(freq_table))) P2 <- arguments[3] else {paste(arguments[3], "not found in header."); q(save="no", status=1)}
      if (is.element(arguments[4], names(freq_table))) P3 <- arguments[4] else {paste(arguments[4], "not found in header."); q(save="no", status=1)}
   } else {
      print("There should be either 1 or 4 arguments after --args")
      q(save="no", status=1)
   }
}

# =====================================================
#  FUNCTIONS
# =====================================================

abba <- function(p1, p2, p3) (1 - p1) * p2 * p3

baba <- function(p1, p2, p3) p1 * (1 - p2) * p3

D.stat <- function(dataframe) (sum(dataframe$ABBA) - sum(dataframe$BABA)) / (sum(dataframe$ABBA) + sum(dataframe$BABA))


# =====================================================
#  MAIN
# =====================================================

ABBA <- abba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
BABA <- baba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
ABBA_BABA_df <- as.data.frame(cbind(ABBA,BABA))
D <- D.stat(ABBA_BABA_df)

write.table(list(D=D), append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

