#!/bin/bash
#
#				2018-07-25
#				----------
#
# On 2018-07-24, Kristyna repeated the filtering of the vcf file. The newly
# filtered vcf should be a common starting point of all downstream analyses.
# Here, I just add the BPF label to the INFO field of the new vcf file. This
# flag is necessary to prepare the input files for the abba/baba test with
# our freq2.awk script.

FILTERED_VCF=/data/kristyna/hedgehog/results_2018/24-07-2018/all.recode.vcf
ADD_FLAG_DIR=../../bin

if [ ! -e flagged.vcf ]; then
   if [ ! -e $FILTERED_VCF ]; then
      echo "Error: Cannot find file $FILTERED_VCF."
      exit
   fi
   if [ ! -e $ADD_FLAG_DIR/add_flag.awk ]; then
      echo "Error: Cannot find script add_flag.awk in $ADD_FLAG_DIR."
      exit
   fi
   gawk -f $ADD_FLAG_DIR/add_flag.awk $FILTERED_VCF > flagged.vcf
fi
