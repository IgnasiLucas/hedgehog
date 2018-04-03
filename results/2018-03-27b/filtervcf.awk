# Usage:
#
#   gawk -f filtervcf.awk [-v MINQ=<INT>] popmap.txt input.vcf > output.vcf
#
# This instance is hard-coded to filter sites where the populations
# specified in the BEGIN block include as many individuals covered as
# indicated there. Change the BEGIN block below to use it in a different
# context. The name of the population map, "popmap.txt" is also hard
# coded. It should have two columns: sample id, and population name.
# The sample names must be exactly the same as in the header of the
# input vcf file. Population names must coincide with those specified
# in the BEGIN block.

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
   MIN["europaeus"]   = 5
   MIN["romanicus"]   = 5
   MIN["concolor"]    = 2
   MIN["Hemiechinus"] = 1

}(FILENAME == "popmap.txt"){
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
