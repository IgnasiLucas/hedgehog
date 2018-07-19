# This script adds one field to the information column of a vcf file. The new
# field is called "Binary presence flag" or BPF, and it is vector of as many
# integers as alleles are present in the variable site. For each allele, the
# bits that compose the integer indicate the presence of that allele in the
# sample. For example, a SNP with alleles A and C genotyped in 6 individuals
# may show allele C in the first one and in the last two. Thus, the bits should
# be: 11001, which in decimal notation is 1+16+32=49. The second allele may be
# present in all individuals, but the last one: 011111 = 62. Note that the "first"
# bit is on the right, and subsequent ones are added to the left, as it should
# be. However, the samples are numbered from left to right in the vcf file.
#
# This is not practical for hundreds of individuals, and impossible for more than
# 1022.
#
# Use with caution, since it involves a modification of a vcf file, and it has
# not been tested extensively.

BEGIN{
   OFS = "\t"
   WRITTEN = 0
}
(/^#/){
   if ((/^##FORMAT/) && (WRITTEN == 0)){
      print "##INFO=<ID=BPF,Number=A,Type=Integer,Description=\"Binary presence flag, indicates what samples each allele is present in.\">"
      WRITTEN = 1
   }
   print
}(/^[^#]/){
   delete FLAG
   FLAG[0] = 0

   # I take the identity of the alleles from the list of alternatives (field 5)
   # and from the reference (field 4). Note they are indexed in the ALLELES
   # array in the same order as in the VCF file.
   split($5, ALLELES, /,/)
   ALLELES[0] = $4

   # GT is the position in the individual information string of the genotype.
   GT = 0
   split($9,ORDER,/:/)
   for (o in ORDER) {
      if (ORDER[o] == "GT") GT = o
   }


   for (i = 0; i <= (NF - 10); i++) {
      split($(i+10), FORMAT, /:/)
      # note that below, 'a' is the index in the array ALLELES: 0, 1...
      # which is how alleles are represented in the GT field: 0/0, 0/1, 1/1...
      # note also that the first individual is the 0th.
      for (a in ALLELES) {
         if (FORMAT[GT] ~ a) FLAG[a] += 2^i
      }
   }

   BPF = ";BPF=" FLAG[0]
   for (a = 1; a < length(ALLELES); a++) {
      BPF = BPF "," FLAG[a]
   }

   $8 = $8 BPF
   print $0
}
