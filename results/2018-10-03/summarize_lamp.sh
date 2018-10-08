#!/bin/bash
#
# Usage: summarize_lamp.sh <dir> <dict>
#
# where <dir> is the individual's directory, that includes the lamp
# results in separate, contig-specific folders. <dict> is the file
# with the contig lengths.
#
# It generates a table with the ancestry percentages, contig length
# and number of snps.

for contig in `find $1 -mindepth 1 -type d -exec basename '{}' \;`; do
   if [ -e $1/$contig/average-anc.txt ] && [ -e $1/$contig/statistics.txt ]; then
      LENGTH=`grep $contig $2 | cut -f 2`
      ANC1=`cut -f 2 $1/$contig/average-anc.txt`
      ANC2=`cut -f 3 $1/$contig/average-anc.txt`
      SNPS=`grep "^SNPs = " $1/$contig/statistics.txt | cut -d " " -f 3`
      echo -e "$contig\t$LENGTH\t$SNPS\t$ANC1\t$ANC2"
   fi
done
