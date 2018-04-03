# Usage:
#
#   gawk [-v OUTGROUP=<str>] [-v P1=<str>] [-v P2=<str>] [-v P3=<str>] -f freq.awk popmap.txt input.vcf
#
# The input vcf must include the binary presence flags.
# The names of the populations and samples must coincide in the
# vcf file and in the population mapping file, which must be
# named 'popmap.txt'. The popmap.txt file may contain more
# populations than the four needed for the abbababa test. In any
# case, the outgroup and populations P1, P2, and P3 must be
# defined either in the command line or in the BEGIN section
# below.

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
   if (OUTGROUP == "") OUTGROUP = "Hemiechinus"
   if (P1 == "") P1 = "romanicus"
   if (P2 == "") P2 = "europaeus"
   if (P3 == "") P3 = "concolor"
   if (MIN1 == "") MIN1 = 3
   if (MIN2 == "") MIN2 = 3
   if (MIN3 == "") MIN3 = 3
   if (MINOUT == "") MINOUT = 1
   POPLIST[P1] = 1
   POPLIST[P2] = 1
   POPLIST[P3] = 1
}(FILENAME == "popmap.txt"){
   POPULATION[$1] = $2
}(/^#CHROM/){
   print "scaffold\tposition\t" P1 "\t" P2 "\t" P3 "\t" OUTGROUP
   for (i = 0; i <= (NF - 10); i++) {
      BINARY_MASK[POPULATION[$(i + 10)]] += 2^i
   }
}((FILENAME != "popmap.txt") && (/^[^#]/)){
   split($8, INFO, /;/)
   for (i in INFO) {
      if (INFO[i] ~ /^BPF=/){
         split(INFO[i], BPF, /=|,/)
      }
   }
   if ((countbits(and(BINARY_MASK[OUTGROUP], BPF[2])) >= MINOUT) && (and(BINARY_MASK[OUTGROUP], BPF[3]) == 0)) {
   # The reference allele is ancestral
      for (pop in POPLIST) {
         AA[pop] = and(and(BINARY_MASK[pop],       BPF[2]), compl(BPF[3]))
         AB[pop] = and(and(BINARY_MASK[pop],       BPF[2]),       BPF[3])
         BB[pop] = and(and(BINARY_MASK[pop], compl(BPF[2])),      BPF[3])
         NAA[pop] = countbits(AA[pop])
         NBB[pop] = countbits(BB[pop])
         NAB[pop] = countbits(AB[pop])
         FREQ[pop] = (2 * NBB[pop] + NAB[pop]) / (2 * (NAA[pop] + NAB[pop] + NBB[pop]))
      }
      if (((NAA[P1] + NAB[P1] + NBB[P1]) >= MIN1) && \
          ((NAA[P2] + NAB[P2] + NBB[P2]) >= MIN2) && \
          ((NAA[P3] + NAB[P3] + NBB[P3]) >= MIN3)) {
         printf "%s\t%u\t%.4f\t%.4f\t%.4f\t%.4f\n", $1, $2, FREQ[P1], FREQ[P2], FREQ[P3], 0
      }

   } else {
      if ((and(BINARY_MASK[OUTGROUP], BPF[2]) == 0) && (countbits(and(BINARY_MASK[OUTGROUP], BPF[3])) >= MINOUT)) {
      # The alternative allele is ancestral
         for (pop in POPLIST) {
            AA[pop] = and(and(BINARY_MASK[pop], compl(BPF[2])),      BPF[3])
            AB[pop] = and(and(BINARY_MASK[pop],       BPF[2]),       BPF[3])
            BB[pop] = and(and(BINARY_MASK[pop],       BPF[2]), compl(BPF[3]))
            NAA[pop] = countbits(AA[pop])
            NAB[pop] = countbits(AB[pop])
            NBB[pop] = countbits(BB[pop])
            FREQ[pop] = (2 * NBB[pop] + NAB[pop]) / (2 * (NAA[pop] + NAB[pop] + NBB[pop]))
         }
         if (((NAA[P1] + NAB[P1] + NBB[P1]) > MIN1) && \
             ((NAA[P2] + NAB[P2] + NBB[P2]) > MIN2) && \
             ((NAA[P3] + NAB[P3] + NBB[P3]) > MIN3)) {
            printf "%s\t%u\t%.4f\t%.4f\t%.4f\t%.4f\n", $1, $2, FREQ[P1], FREQ[P2], FREQ[P3], 0
         }
      }
   }
}
