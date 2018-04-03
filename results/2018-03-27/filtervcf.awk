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
   if (REQUIRED == 0) {
      print "You have to set the variable REQUIRED."
      print "Only sites including all samples specified in the REQUIRED variable"
      print "will be output."
      exit
   }
}(/^#/){
   print $0
}(/^[^#]/){
   split($8, INFO, /;/)
   for (i in INFO) {
      if (INFO[i] ~ /^BPF=/){
         split(INFO[i], BPF, /=|,/)
      }
   }
   WITHDATA = or(BPF[2], BPF[3]) # dismissing 3rd and subsequent alleles
#   FILTERED = countbits(and(REQUIRED, WITHDATA)) # this counts how many of the required samples have either the first or the second allele
   if (REQUIRED == and(REQUIRED, WITHDATA)){
      print $0
   }
}
