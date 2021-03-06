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
   MINEURO   = 12
   MINROMA   = 12
   MINCONCO  = 2
   EUROBITS  = 251293811610192
   ROMABITS  = 566536679968143
   CONCOBITS = 299273322233856
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
   EUROPAEUS = countbits(and(EUROBITS, WITHDATA))
   ROMANICUS = countbits(and(ROMABITS, WITHDATA))
   CONCOLOR  = countbits(and(CONCOBITS, WITHDATA))
   if ((EUROPAEUS >= MINEURO) && (ROMANICUS >= MINROMA) && (CONCOLOR >= MINCONCO)) {
      print $0
   }
}
