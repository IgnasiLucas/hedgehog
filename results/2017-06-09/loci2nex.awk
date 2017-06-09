(/^\/\//){
   PIS   = gsub(/*/,"*")
   split($NF, A, /\|/)
   OUTFILE = sprintf("loc%s.nex", A[2])
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
         if (NAME ~ /_1$/) {
            if (SEQ[NAME] != SEQ[substr(NAME,1,length(NAME)-2) "_0"]) {
               printf "    %-*s  %s\n", LONGEST, NAME, SEQ[NAME] >OUTFILE
            }
         } else {
            printf "    %-*s  %s\n", LONGEST, NAME, SEQ[NAME] >OUTFILE
         }
      }
      print "  ;" >OUTFILE
      print "end;" >OUTFILE
   }
   close(OUTFILE)
   delete(SEQ)
   NTAX = 0
}(/^[^\/]/){
   SEQ[$1] = $2
   NTAX++
   NCHAR = length($2)
}
