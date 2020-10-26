# To turn decimal numbers to binary
#
#    awk -f dec2bin.awk <file>
#
# With decimal number in first field per line in input file.
# This was taken from an answer by James Brown in StackOverflow.

function d2b(d,  b) {
   while(d) {
      b = d%2 b
      d = int(d/2)
   }
   return(b)
}
function d2b2(d,z,  b){
   b = ""
   while(z + 1){
      if (d >= 2^z) {
         b = b 1
         d = d - 2^z
         z = z-1
      } else {
         b = b 0
         z = z-1
      }
   }
   return(b)
}

BEGIN{
   if (z == 0) z = 76
}{
#   S = sprintf("%*s", z, d2b($1))
#   gsub(/ /,"0",S)
#   print S
   print d2b2($1,z-1)
}
