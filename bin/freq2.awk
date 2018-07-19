# Usage:
#
#   gawk [-v <VAR>=<VALUE> ...] -f freq.awk popmap.txt input.vcf
#
# The purpose of this script is to parse a vcf file with the BPF information
# tag and print the derived allele frequencies of 4 or 5  populations at variable
# sites where the last population (the outgroup) is fixed for the ancestral
# allele. It will also print for each site how many individuals from each
# population were available to estimate allele frequencies, if requested.
#
# Options are passed as variable definitions in the command line, with the
# "-v" argument of awk. The available options are:
#
#  -v P1=<str>		Name of population 1. Default: "romanicus".
#  -v P2=<str>		Name of population 2. Default: "concolor".
#  -v P3=<str>		Name of population 3. Default: "europaeus".
#  -v P4=<str>		Name of population 4. Default: none.
#  -v OUTGROUP=<str>	Name of outgroup population. Def.: "Hemiechinus".
#  -v PRINT_N=<int>	Request printing the number of individuals used if <int> > 0.
#  -v MIN1=<int>	Minimum number of individuals in population 1. Def.: 5.
#  -v MIN2=<int>	Minimum number of individuals in population 2. Def.: 5.
#  -v MIN3=<int>	Minimum number of individuals in population 3. Def.: 5.
#  -v MINOUT=<int>	Minimum number of individuals in outgroup. Default: 1.
#
# The input vcf must include the binary presence flags.
# The names of the populations and samples must coincide in the
# vcf file and in the population mapping file, which must be
# named 'popmap.txt'. The popmap.txt file may contain more
# populations than the four needed for the abbababa test. In any
# case, the outgroup and populations P1, P2, and P3 must be
# defined either in the command line or in the BEGIN section
# below.

function countbits(bits,    data){
   if (bits == 0) return 0
   for (; bits != 0; bits = rshift(bits, 1)) data += and(bits, 1)
   return data
}

BEGIN{
   if (OUTGROUP == "") OUTGROUP = "Hemiechinus"
   if (P1 == "") P1 = "romanicus"
   if (P2 == "") P2 = "concolor"
   if (P3 == "") P3 = "europaeus"
   if (MIN1 == "") MIN1 = 5
   if (MIN2 == "") MIN2 = 5
   if (MIN3 == "") MIN3 = 5
   if (MIN4 == "") MIN4 = 5
   if (MINOUT == "") MINOUT = 1
   POPLIST[1] = P1
   POPLIST[2] = P2
   POPLIST[3] = P3
   if (P4 != "") {
      POPLIST[4] = P4
   }
}(FILENAME ~ /popmap/){
   POPULATION[$1] = $2
   POP_IN_MAP[$2]++
   IND_IN_MAP[$1]++
}((FILENAME ~ /vcf$/) && (/^#CHROM/)){
   for (i = 0; i <= (NF - 10); i++) {
      if ($(i + 10) in IND_IN_MAP) {
         BINARY_MASK[POPULATION[$(i + 10)]] += 2^i
      } else {
         print "Individual " $(i + 10) " in " FILENAME " from unknown population."
         exit 1
      }
   }
   for (i in POPLIST) {
      if (!(POPLIST[i] in POP_IN_MAP)) {
         print "Population " POPLIST[i] " is not found in the population map file."
         exit 1
      }
   }
   if (4 in POPLIST) {
      line = sprintf("scaffold\tposition\t%s\t%s\t%s\t%s\t%s", P1, P2, P3, P4, OUTGROUP)
   } else {
      line = sprintf("scaffold\tposition\t%s\t%s\t%s\t%s", P1, P2, P3, OUTGROUP)
   }
   if (PRINT_N > 0) {
      for (i = 1; i <= length(POPLIST); i++) {
         line = line "\tN" i
      }
   }
   print line
   line = ""
}((FILENAME ~ /vcf$/) && (/^[^#]/)){
   split($8, INFO, /;/)
   for (i in INFO) {
      if (INFO[i] ~ /^BPF=/){
         split(INFO[i], BPF, /=|,/)
      }
   }
   if ((countbits(and(BINARY_MASK[OUTGROUP], BPF[2])) >= MINOUT) && (and(BINARY_MASK[OUTGROUP], BPF[3]) == 0)) {
   # The reference allele is ancestral
      for (pop in POPLIST) {
         AA[pop] = and(and(BINARY_MASK[POPLIST[pop]],       BPF[2]), compl(BPF[3]))
         AB[pop] = and(and(BINARY_MASK[POPLIST[pop]],       BPF[2]),       BPF[3])
         BB[pop] = and(and(BINARY_MASK[POPLIST[pop]], compl(BPF[2])),      BPF[3])
         NAA[pop] = countbits(AA[pop])
         NBB[pop] = countbits(BB[pop])
         NAB[pop] = countbits(AB[pop])
         N[pop] = NAA[pop] + NAB[pop] + NBB[pop]
         if (N[pop] > 0) {
            FREQ[pop] = (2 * NBB[pop] + NAB[pop]) / (2 * N[pop])
         } else {
            FREQ[pop] = 999.9999
         }
      }
   } else {
      if ((and(BINARY_MASK[OUTGROUP], BPF[2]) == 0) && (countbits(and(BINARY_MASK[OUTGROUP], BPF[3])) >= MINOUT)) {
      # The alternative allele is ancestral
         for (pop in POPLIST) {
            AA[pop] = and(and(BINARY_MASK[POPLIST[pop]], compl(BPF[2])),      BPF[3])
            AB[pop] = and(and(BINARY_MASK[POPLIST[pop]],       BPF[2]),       BPF[3])
            BB[pop] = and(and(BINARY_MASK[POPLIST[pop]],       BPF[2]), compl(BPF[3]))
            NAA[pop] = countbits(AA[pop])
            NAB[pop] = countbits(AB[pop])
            NBB[pop] = countbits(BB[pop])
            N[pop] = NAA[pop] + NAB[pop] + NBB[pop]
            if (N[pop] > 0) {
               FREQ[pop] = (2 * NBB[pop] + NAB[pop]) / (2 * (NAA[pop] + NAB[pop] + NBB[pop]))
            } else {
               FREQ[pop] = 999.9999
            }
         }
      } else {
         next
      }
   }
   if ((N[1] >= MIN1) && (N[2] >= MIN2) && (N[3] >= MIN3)) {
      line = sprintf("%s\t%u\t%.4f\t%.4f\t%.4f", $1, $2, FREQ[1], FREQ[2], FREQ[3])
   } else {
      next
   }
   if (4 in POPLIST) {
      line = line sprintf("\t%.4f\t%.4f", FREQ[4], 0)
   } else {
      line = line sprintf("\t%.4f", 0)
   }
   if (PRINT_N > 0) {
      for (i = 1; i <= length(POPLIST); i++) {
         line = line "\t" N[i]
      }
   }
   print line
   line = ""
}
