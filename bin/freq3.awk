# Usage:
#
#   gawk --bignum [-v <VAR>=<VALUE> ...] -f freq.awk popmap.txt input.vcf
#
# The purpose of this script is to parse a vcf file with the BPF information
# tag and print the derived allele frequencies of X  populations at variable
# sites where the last population (the outgroup) is fixed for the ancestral
# allele. It will also print for each site how many individuals from each
# population were available to estimate allele frequencies, if requested.
#
# Options are passed as variable definitions in the command line, with the
# "-v" argument of awk. The available options are:
#
#  -v POPS=<pop1,pop2,pop3...>	Names of populations 1, 2, 3... Comma-
#                               delimited. Last population assumed to be
#                               the outgroup. Default: "".
#  -v MINS=<1,5,5...> 		Minimum number of individuals per population,
#                               in the same order as the populations.
#  -v PRINT_N=<int>	Request printing the number of individuals used if <int> > 0.
#
# The input vcf must include the binary presence flags.
# The names of the populations and samples must coincide in the
# vcf file and in the population mapping file, which must be
# named 'popmap.txt'. The popmap.txt file may contain more
# populations than the ones requested.

function countbits(bits,    data){
   if (bits == 0) return 0
   for (; bits != 0; bits = rshift(bits, 1)) data += and(bits, 1)
   return data
}

BEGIN{
   if (POPS == "") {
      print "Populations undefined."
      exit 1
   }
   if (MINS == "") {
      print "Minimum numbers of individuals per population undefined."
      exit 1
   }
   split(POPS,POPLIST,",")
   split(MINS,MINLIST,",")
   NUMPOPS = length(POPLIST)
   OUTGROUP = POPLIST[NUMPOPS]
   MINOUT = MINLIST[NUMPOPS]
   if (length(MINLIST) != NUMPOPS) {
      print "Different number of populations than of minimum numbers of individuals."
      exit 1
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
   line = "scaffold\tposition"
   for (pop = 1; pop <= NUMPOPS; pop++){
      line = line "\t" POPLIST[pop]
   }
   if (PRINT_N > 0) {
      for (i = 1; i <= NUMPOPS; i++) {
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
   GOOD = 0
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
         if (N[pop] >= MINLIST[pop]) {
            FREQ[pop] = (2 * NBB[pop] + NAB[pop]) / (2 * N[pop])
            GOOD += 1
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
            if (N[pop] >= MINLIST[pop]) {
               FREQ[pop] = (2 * NBB[pop] + NAB[pop]) / (2 * (NAA[pop] + NAB[pop] + NBB[pop]))
               GOOD += 1
            }
         }
      } else {
         next
      }
   }
   if (GOOD == NUMPOPS){
      line = $1 "\t" $2
      for (pop =1; pop <= NUMPOPS; pop++) {
         line = line "\t" FREQ[pop]
      }
   } else {
      next
   }
   if (PRINT_N > 0) {
      for (pop = 1; pop <= NUMPOPS; pop++) {
         line = line "\t" N[pop]
      }
   }
   print line
   line = ""
}
