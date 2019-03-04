# Usage:
#
#   gawk -f filtervcf.awk [-v MINQ=<INT>] [-v EUR=<int>] [-v ROU=<int>] [-v CON=<int>] [-v HEM=<int>] popmap.txt input.vcf > output.vcf
#
# This simple script filters a vcf file with individuals from four
# different populations, the names of which are hardcoded for the sake of
# the present project: europaeus, roumanicus, concolor, Hemiechinus. The
# main filter selects sites where a minimum number of individuals from each
# population have genotype data (genotype is not "./."). You can specify the
# minimum number of individuals from each population using the options:
#    -v EUR=<int>
#    -v ROU=<int>
#    -v CON=<int>
#    -v HEM=<int>
# Defaults are: 10 europaeus, 10 roumanicus, 4 concolor and 0 Hemiechinus.
# Note that the name of the file specifying the distribution of samples among
# populations must include the string "popmap", and must be given before the
# input vcf file in the command line. The popmap file must have one individual
# per line, with two fields: its name (as used in the vcf file) and the name of
# the population where it belongs (as used in this script).
#
# The script also checks that the quality of the variable site is above a minimum,
# if MINQ is set, via "-v MINQ=<int>". And it also requires sites to be bi-allelic.
#
# The input vcf file must include the binary presence flag (BPF field in the
# INFO section. See add_flag.awk for instructions.

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
   if (EUR == 0) {MIN["europaeus"]   = 0} else {MIN["europaeus"]   = EUR}
   if (ROU == 0) {MIN["roumanicus"]  = 0} else {MIN["roumanicus"]  = ROU}
   if (CON == 0) {MIN["concolor"]    = 0} else {MIN["concolor"]    = CON}
   if (HEM == 0) {MIN["Hemiechinus"] = 0} else {MIN["Hemiechinus"] = HEM}

}(FILENAME ~ /popmap/){
   POPULATION[$1] = $2
}
(/^##/){
   print $0
}(/^#CHROM/){
   print $0
   for (i = 0; i <= (NF - 10); i++) {
      BINARY_MASK[POPULATION[$(i + 10)]] += 2^i
   }
}(/^[^#]/){
   split($8, INFO, /;/)
   for (i in INFO) {
      if (INFO[i] ~ /^BPF=/){
         NUM_ALLELES = split(INFO[i], BPF, /=|,/) - 1
      }
   }
   WITHDATA = or(BPF[2], BPF[3]) # dismissing 3rd and subsequent alleles
   TESTS = 0
   for (population in MIN) {
      TESTS += (countbits(and(BINARY_MASK[population], WITHDATA)) >= MIN[population])
   }
   TESTS += (NUM_ALLELES == 2)
   TESTS += ($6 >= MINQ)
   if (TESTS == length(MIN) + 2) {
      print $0
   }
}
