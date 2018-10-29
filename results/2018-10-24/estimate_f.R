# I separate the estimation of f, the proportion of the genome that is admixed, from the
# estimation of D, in script abba_baba.R. The reason is that the input files may be different,
# and it was getting difficult to add options to the abba_baba.R script.
#
# =====================================================
#  USAGE
# =====================================================
#
# R [options] <estimate_f.R [>outfile] --args <freq_table> [<P1> <P2> <P3a> <P3b> <P3>] 
#
#    options: --save, --no-save or --vanilla. One required.
#
#    freq_table  Tab-separated values of derived allele frequencies in populations P1, P2, P3a
#                and P3b. It must have a header with the population names. Populations P1, P2, P3a
#                and P3b may be asumed to correspond to columns 3, 4, 5 and 6, after 'scaffold' and
#                'position', unless population names are given in a different order in the command
#                line.
#
#    P1 P2 P3a P3b P3
#                If the population names are given as additional arguments, the order in which
#                they are specified will define which columns in freq_table correspond to P1, P2,
#                P3a and P3b, as long as they coincide with the names in the header of freq_table.
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
   if (dim(freq_table)[2] >= 6) {
      P1  <- names(freq_table)[3]
      P2  <- names(freq_table)[4]
      P3a <- names(freq_table)[5]
      P3b <- names(freq_table)[6]
      P3  <- names(freq_table)[7]
   } else {
      # For estimation of f we don't need the scaffold and position.
      if (dim(freq_table)[2] == 4) {
         P1  <- names(freq_table)[1]
         P2  <- names(freq_table)[2]
         P3a <- names(freq_table)[3]
         P3b <- names(freq_table)[4]
         P3  <- names(freq_table)[5]
      } else {
         print("Bad format of freq_table.")
         q(save="no", status=1)
      }
   }
} else {
   if (length(arguments) == 6){
      freq_table <- read.table(arguments[1], header=TRUE, as.is=TRUE)
      if (is.element(arguments[2], names(freq_table))) P1  <- arguments[2] else {paste(arguments[2], "not found in header."); q(save="no", status=1)}
      if (is.element(arguments[3], names(freq_table))) P2  <- arguments[3] else {paste(arguments[3], "not found in header."); q(save="no", status=1)}
      if (is.element(arguments[4], names(freq_table))) P3a <- arguments[4] else {paste(arguments[4], "not found in header."); q(save="no", status=1)}
      if (is.element(arguments[5], names(freq_table))) P3b <- arguments[5] else {paste(arguments[5], "not found in header."); q(save="no", status=1)}
      if (is.element(arguments[6], names(freq_table))) P3  <- arguments[6] else {paste(arguments[6], "not found in header."); q(save="no", status=1)}
   } else {
      print("There should be either 1 or 6 arguments after --args")
      q(save="no", status=1)
   }
}

# =====================================================
#  FUNCTIONS
# =====================================================

abba <- function(p1, p2, p3) (1 - p1) * p2 * p3

baba <- function(p1, p2, p3) p1 * (1 - p2) * p3

abba_d <- function(p1, p2, p3){
   abba(p1, pmax(p2, p3), pmax(p2, p3))
}

baba_d <- function(p1, p2, p3){
   baba(p1, pmax(p2, p3), pmax(p2, p3))
}


# =====================================================
#  MAIN
# =====================================================

ABBA_1_2_3a  <-   abba(freq_table[,P1], freq_table[,P2],  freq_table[,P3a])
BABA_1_2_3a  <-   baba(freq_table[,P1], freq_table[,P2],  freq_table[,P3a])

ABBA_1_3b_3a <-   abba(freq_table[,P1], freq_table[,P3b], freq_table[,P3a])
BABA_1_3b_3a <-   baba(freq_table[,P1], freq_table[,P3b], freq_table[,P3a])

ABBA_1_2_3   <-   abba(freq_table[,P1], freq_table[,P2],  freq_table[,P3])
BABA_1_2_3   <-   baba(freq_table[,P1], freq_table[,P2],  freq_table[,P3])

ABBA_1_D_D   <- abba_d(freq_table[,P1], freq_table[,P2],  freq_table[,P3])
BABA_1_D_D   <- baba_d(freq_table[,P1], freq_table[,P2],  freq_table[,P3])

ABBA_1_3_3   <-   abba(freq_table[,P1], freq_table[,P3],  freq_table[,P3])
BABA_1_3_3   <-   baba(freq_table[,P1], freq_table[,P3],  freq_table[,P3])

f_hom <- (sum(ABBA_1_2_3)   - sum(BABA_1_2_3)) /
         (sum(ABBA_1_3_3)   - sum(BABA_1_3_3))

f     <- (sum(ABBA_1_2_3a)  - sum(BABA_1_2_3a)) /
         (sum(ABBA_1_3b_3a) - sum(BABA_1_3b_3a))

f_d   <- (sum(ABBA_1_2_3)   - sum(BABA_1_2_3)) /
         (sum(ABBA_1_D_D)   - sum(BABA_1_D_D))

write.table(list(P1=P1, P2=P2, P3a=P3a, P3b=P3b, P3=P3, f_hom=f_hom, f_d=f_d, f=f), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)


