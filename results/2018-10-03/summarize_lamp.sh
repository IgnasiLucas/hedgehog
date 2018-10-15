#!/bin/bash
#
# Usage: summarize_lamp.sh <dir> <gen> <dict>
#
# where <dir> is the individual's directory, that includes the lamp
# results in separate, contig-specific folders. <dict> is the file
# with the contig lengths.
#
# It generates a table with the ancestry percentages, contig length
# and number of snps.

for contig in `find $1 -mindepth 1 -name 'NW*' -type d -exec basename '{}' \;`; do
   if [ -e $(printf "%s/%s/gen%06i/average-anc.txt" $1 $contig $2) ] && [ -e $(printf "%s/%s/gen%06i/statistics.txt" $1 $contig $2) ]; then
      LENGTH=`grep $contig $3 | cut -f 2`
      ANC1=`cut -f 2 $(printf "%s/%s/gen%06i/average-anc.txt" $1 $contig $2)`
      ANC2=`cut -f 3 $(printf "%s/%s/gen%06i/average-anc.txt" $1 $contig $2)`
      SNPS=`grep "^SNPs = " $(printf "%s/%s/gen%06i/statistics.txt" $1 $contig $2) | cut -d " " -f 3`
      echo -e "$contig\t$LENGTH\t$SNPS\t$ANC1\t$ANC2"
   fi
done
