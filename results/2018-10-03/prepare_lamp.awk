# This awk script taks as input files first the population mapping
# file (usually, popmap.txt) and a vcf file with binary presence flags.
# The individuals that must be analysed must be passed as a comma-
# separated list of sample names to variable IND. Similarly, the ancestral
# population 1 and 2 of each admixed individual must be specified as
# comma-separated files in variables POP1 and POP2. And their proportions
# of admixture must be specified in Q1 and Q2, also as comma-separated
# values. For example:
#
#   gawk -v -v IND="Er37_SK27,Er55_AU7,Er26_JUG4,Er27_SK32" \
#        -v POP1="roumanicus,roumanicus,roumanicus,roumanicus" \
#        -v POP2="europaeus,europaeus,concolor,concolor" \
#        -v Q1="0.747088,0.951943,0.960462,0.976133" \
#        -v Q2="0.252902,0.048047,0.039528,0.023857" \
#        -f prepare_lamp.awk popmap.txt <VCF>
#
# The outputs will be written in <individual>/<contig>/ folders, and they
# include: a genotypes file for the individual in all variable positions
# of the contig (geno.txt); two allele frequencies files, one for each
# ancestral population (Anc1.P.txt, Anc2.P.txt), and a positions file
# (positions.txt).

function bits2str(bits, max,        data, mask){
   if (bits == 0)
      return "0"

   mask = 1
   for (; bits != 0; bits = rshift(bits, 1))
      data = (and(bits, mask) ? "1" : "0") data

   while ((length(data) % max) != 0)
      data = "0" data

   return data
}

function countbits(bits,    data){
   if (bits == 0) return 0
   for (; bits != 0; bits = rshift(bits, 1)) data += and(bits, 1)
   return data
}

BEGIN{
   N1 = split(IND,INDA,/,/)
   for (ind in INDA) {
      # The INDB array has sample names as keys.
      INDB[INDA[ind]] = ind
      SUM_SIZE_1[INDA[ind]] = 0
      SUM_SIZE_2[INDA[ind]] = 0
      NUM_SITES[INDA[ind]] = 0
   }
   N2 = split(POP1,POP1A,/,/)
   N3 = split(POP2,POP2A,/,/)
   if ((N1 != N2) || (N2 != N3)){
      print "arrays are not the same size"
      exit 1
   }
}(FILENAME ~ /popmap/){
   POPULATION[$1] = $2
   ALL_POPS[$2] = 1
}(FILENAME ~ /\.vcf$/){
   if ($1 ~ /^#CHROM/){
      print $0
      for (i = 0; i <= (NF - 10); i++) {
         BINARY_MASK[POPULATION[$(i + 10)]] += 2^i
         BINARY_IND[$(i + 10)] = 2^i
      }
      for (pop in BINARY_MASK) {printf("% 10s\t% 13i\t%s\t% 3i\n", pop, BINARY_MASK[pop], bits2str(BINARY_MASK[pop], 50), countbits(BINARY_MASK[pop]))}
   }
   if ($0 ~ /^[^#]/){
      split($8, INFO, /;/)
      for (i in INFO) {
         if (INFO[i] ~ /^BPF=/){
            split(INFO[i], BPF, /=|,/)
         }
      }
      CONTIG = $1
      POS    = $2
      WITHDATA =  or(BPF[2], BPF[3])
      HETERO   = and(BPF[2], BPF[3])
      HOMO_REF = and(BPF[2], compl(BPF[3]))
      HOMO_ALT = and(compl(BPF[2]), BPF[3])
      for (ind in INDA) {
         # ind is the individual's index in INDA. It's name is INDA[ind]
         # Genotype is the number of alternative alleles (BPF[3])
         if (and(WITHDATA, BINARY_IND[INDA[ind]])) {
#
#                                    1 if ALT is present,                           1 if REF is absent,
#                                    0 if ALT is absent                             0 if REF is present
#                       ______________________^___________________     ______________________^____________________
#                      /                                          \   /                                           \
            GENOTYPE = countbits(and(BPF[3], BINARY_IND[INDA[ind]])) + countbits(and(compl(BPF[2]), BINARY_IND[INDA[ind]]))
         } else {
            GENOTYPE = -1
         }
         # Below, I exclude the focus individual from its alleged population
         # before estimating allele frequencies.
         HOMO_ALT_1 = countbits(and(HOMO_ALT, and(BINARY_MASK[POP1A[ind]], compl(BINARY_IND[INDA[ind]]))))
         HETERO_1   = countbits(and(HETERO,   and(BINARY_MASK[POP1A[ind]], compl(BINARY_IND[INDA[ind]]))))
         TOTAL_1    = countbits(and(WITHDATA, and(BINARY_MASK[POP1A[ind]], compl(BINARY_IND[INDA[ind]]))))
         if (TOTAL_1 > 0) {
            FREQ_1 = (HETERO_1 + 2 * HOMO_ALT_1) / (2 * TOTAL_1)
            SUM_SIZE_1[INDA[ind]] += TOTAL_1
         } else {
            FREQ_1 = -1
         }
         HOMO_ALT_2 = countbits(and(HOMO_ALT, and(BINARY_MASK[POP2A[ind]], compl(BINARY_IND[INDA[ind]]))))
         HETERO_2   = countbits(and(HETERO,   and(BINARY_MASK[POP2A[ind]], compl(BINARY_IND[INDA[ind]]))))
         TOTAL_2    = countbits(and(WITHDATA, and(BINARY_MASK[POP2A[ind]], compl(BINARY_IND[INDA[ind]]))))
         if (TOTAL_2 > 0) {
            FREQ_2 = (HETERO_2 + 2 * HOMO_ALT_2) / (2 * TOTAL_2)
            SUM_SIZE_2[INDA[ind]] += TOTAL_2
         } else {
            FREQ_2 = -1
         }
         if ((FREQ_1 > -1) && (FREQ_2 > -1)) {
            print POS      >INDA[ind] "/" CONTIG "/positions.txt"
            printf("%i ", GENOTYPE)  >INDA[ind] "/" CONTIG "/geno.txt"
            printf("%.4f\n", FREQ_1) >INDA[ind] "/" CONTIG "/Anc1.P.txt"
            printf("%.4f\n", FREQ_2) >INDA[ind] "/" CONTIG "/Anc2.P.txt"
            NUM_SITES[INDA[ind]] += 1
         }
      }
   }
}END{
   print "\n\nIndividual\tAncestral_1\tAncestral_2\tNum.Sites\tSmpl.Size_1\tSmpl.Size_2"
   for (ind in INDA) {
      printf("% 9s\t% 10s\t% 10s\t% 9i\t% 11.2f\t% 11.2f\n", INDA[ind], POP1A[ind], POP2A[ind], NUM_SITES[INDA[ind]], SUM_SIZE_1[INDA[ind]] / NUM_SITES[INDA[ind]], SUM_SIZE_2[INDA[ind]] / NUM_SITES[INDA[ind]])
   }
}
