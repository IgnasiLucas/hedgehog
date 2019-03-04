# Usage:
#
#   gawk [-v <VAR>=<VALUE> ...] -f vcf2MK.awk popmap.txt input.vcf
#
# The purpose of this script is to parse a vcf file with the BPF information
# tag and print the input file for an asymptotic McDonald-Kretiman test:
# http://benhaller.com/messerlab/asymptoticMK.html.
#
# First, it reads a population file, the name of which must include the
# string 'popmap'. The popmap file is a space- or tab-delimited file with
# two columns: the sample name (the same as in the input vcf file) and the
# name of the population to which the sample belongs. Then, it classifies variants
# as synonymous and non-synonymous (otherwise discarded), and as either
# variable in a focal population or fixed with respect to a sister population.
# It relies on the presence of the flag ANN in the INFO field.
#
# Options are passed as variable definitions in the command line, with the
# "-v" argument of awk. The available options are:
#
#  -v FOCAL=<str>       Name of the focal populaton. Default: "roumanicus".
#  -v SISTER=<str>      Name of the sister population. Default: "europaeus".
#  -v OUTPUT1=<str>     Name of the file to write polymorphism data. Default: "z1.txt".
#  -v OUTPUT2=<str>     Name of the file to write divergence data. DEfault: "z2.txt".

function countbits(bits,    data){
   if (bits == 0) return 0
   for (; bits != 0; bits = rshift(bits, 1)) data += and(bits, 1)
   return data
}

function altfreq(ref, alt, pop,     N, P, H, Q){
   N = countbits(and(pop, or(ref, alt))) # Number of individuals with either ref or alt alleles.
   P = countbits(and(pop, and(ref, compl(alt)))) # Number of individuals homozygous for the reference.
   H = countbits(and(pop, and(ref, alt))) # Number of heterozygous individuals.
   Q = countbits(and(pop, and(alt, compl(ref)))) # Number of individuals homozygous for the alternative.
   if (P + H + Q != N) {
      print "The fucntion altfreq does not work. Sorry."
      exit 1
   }
   return (2 * Q + H) / (2 * N)
}

BEGIN{
   if (FOCAL == "") FOCAL = "roumanicus"
   if (SISTER == "") SISTER = "europaeus"
   if (OUTPUT1 == "") OUTPUT1 = "z1.txt"
   if (OUTPUT2 == "") OUTPUT2 = "z2.txt"
   DIVERGENCE["synonymous"] = 0.0
   DIVERGENCE["missense"] = 0.0
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
   if (!(FOCAL in POP_IN_MAP)) {
      print "Population " FOCAL " is not present in the population map file."
      exit 1
   }
   if (!(SISTER in POP_IN_MAP)) {
      print "Population " SISTER " is not present in the population map file."
      exit 1
   }
}((FILENAME ~ /vcf$/) && (/^[^#]/)){
   EFFECT = ""
   if ($8 ~ /synonymous_variant/) {
      EFFECT = "synonymous"
   }
   if ($8 ~ /missense_variant/) {
      EFFECT = "missense"
   }
   if (EFFECT == "") next
   split($8, INFO, /;/)
   for (i in INFO) {
      if (INFO[i] ~ /^BPF=/){
         split(INFO[i], BPF, /=|,/)   # BPF[2] is the presence of the reference allele, and BPF[3], that of the alternative.
      }
   }
   ALT_FOCAL_FREQ = altfreq(BPF[2], BPF[3], BINARY_MASK[FOCAL])
   ALT_SISTER_FREQ = altfreq(BPF[2], BPF[3], BINARY_MASK[SISTER])
   if (ALT_FOCAL_FREQ + 0.0 == 0.0) {
      DIVERGENCE[EFFECT] += ALT_SISTER_FREQ
   } else {
      if (ALT_FOCAL_FREQ + 0.0 < 1.0) {
         POLYMORPHISM[sprintf("%.1f", ALT_FOCAL_FREQ), EFFECT] += 1
      } else {
         DIVERGENCE[EFFECT] += 1.0 - ALT_SISTER_FREQ
      }
   }
}END{
   for (x = 0.0; x <= 1.0; x += 0.1) {
      print x "\t" POLYMORPHISM[sprintf("%.1f",x), "missense"] + 0 "\t" POLYMORPHISM[x, "synonymous"] + 0 >> OUTPUT1
   }
   print "Synonymous fixed differences: " DIVERGENCE["synonymous"] > OUTPUT2
   print "Non-synonymous fixed differences: " DIVERGENCE["missense"] >> OUTPUT2
}
