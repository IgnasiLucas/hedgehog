#!/bin/bash
#
# Below I use only primary alignments of second reads to
# count as duplicates those that map to the exact same
# position. Then, I print a table of duplicity frequencies.
#
# 260 = 256 + 4: Secondary alignments and unmapped.
# 128: Read 2

samtools view -h -u -F 260 -f 128 $1 2> $2'_'$3.view.err | \
samtools sort -O SAM -T $2 - 2> $2'_'$3.sort.err | \
gawk 'BEGIN{
   CONTIG = "NW_006803924.1"
   POSITION = 1
   NUM_DUPS = 0
   MAX_DUPS = 1
}(/^[^@]/){
   if (($3 == CONTIG) && ($4 == POSITION)) {
      NUM_DUPS++
   } else {
      if (NUM_DUPS > MAX_DUPS) MAX_DUPS = NUM_DUPS
      F[NUM_DUPS]++
      NUM_DUPS = 1
      CONTIG = $3
      POSITION = $4
   }
}END{
   if (NUM_DUPS > MAX_DUPS) MAX_DUPS = NUM_DUPS
   F[NUM_DUPS]++
   SUM = 0
   for (n = 1; n <= MAX_DUPS; n++) {
      SUM += n * F[n]
      print n "\t" F[n] "\t" SUM
   }
}'

