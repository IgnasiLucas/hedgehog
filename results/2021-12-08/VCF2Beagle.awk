(/^#CHROM/){
   TRANSLATE["A"] = 0
   TRANSLATE["C"] = 1
   TRANSLATE["G"] = 2
   TRANSLATE["T"] = 3
   GLfield = 0
   HEADER = "marker\tallele1\tallele2\t"
   for (i = 10; i <= NF; i++) {
      HEADER = HEADER "\t" $i "\t" $i "\t" $i
   }
   print HEADER
}((/^[^#]/) && ($5 !~ /,/) && ($4 ~ /^[ACGT]$/)){
   MARKER  = $1 "_" $2
   ALLELE1 = TRANSLATE[$4]
   ALLELE2 = TRANSLATE[$5]
   LINE = MARKER "\t" ALLELE1 "\t" ALLELE2
   if (GLfield == 0) {
      split($9, FORMAT, /:/)
      for (field in FORMAT) {
         if (FORMAT[field] == "GL") GLfield = field
      }
   }
   for (i = 10; i <= NF; i++) {
      split($i, SAMPLEINFO, /:/)
      split(SAMPLEINFO[GLfield], GL, /,/)
      LINE = LINE sprintf("\t%.6f\t%.6f\t%.6f", 10^GL[1], 10^GL[2], 10^GL[3])
   }
   print LINE
}
