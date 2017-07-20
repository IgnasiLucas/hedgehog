(/^\/\//){
   PIS   = gsub(/*/,"*")
   split($NF, A, /\|/)
   OUTFILE = sprintf("loc%s.nex", A[2])
   # I select loci with more than 2 parsimony informative sites.
   if (PIS > 2) {
      print "#nexus" >OUTFILE
      print "begin data;" >OUTFILE
      printf "  dimensions ntax=%u nchar=%u;\n", NTAX, NCHAR >OUTFILE
      print "  format datatype=dna missing=N gap=-;" >OUTFILE
      print "  matrix" >OUTFILE
      LONGEST = 0
      for (NAME in SEQ) {
         if (length(NAME) > LONGEST) LONGEST = length(NAME)
      }
      for (NAME in SEQ) {
         printf "    %-*s  %s\n", LONGEST, NAME, SEQ[NAME] >OUTFILE
      }
      print "  ;" >OUTFILE
      print "end;" >OUTFILE
   }
   close(OUTFILE)
   delete(SEQ)
   NTAX = 0
}(/^[^\/]/){
   # Note that below I exclude the hybrid individual, labeled as  "Er37".
   if ($1 ~ /_1$/) {
      if ($1 !~ /Er37/) {
         if ($2 != SEQ[substr($1, 1, length($1)-2) "_0"]) {
            SEQ[$1] = $2
            NTAX++
            NCHAR = length($2)
         }
      }
   } else {
      if ($1 !~ /Er37/) {
         SEQ[$1] = $2
         NTAX++
         NCHAR = length($2)
      }
   }
}
