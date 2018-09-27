#!/bin/bash
#
#				2018-09-26
#				==========
#
# I want to re-do the filtering of the vcf file. We need to insert the binary
# presence flag, and then make sure that we select for analysis sites where at least
# four individuals from E. concolor have data.

GZVCF=/data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz
ADD_FLAG_DIR=../../bin
INDLIST=/data/kristyna/hedgehog/results_2018/24-07-2018/popmap

if [ ! -e IndList.txt ]; then
   # This file is just the list of samples, after removing 5 with low quality data. 45 left.
   cp $INDLIST IndList.txt
fi
if [ ! -e r10e10c5.vcf ]; then
   if [ ! -e flagged.vcf ]; then
      vcftools --gzvcf $GZVCF \
            --stdout \
            --keep IndList.txt \
            --maf 0.01 \
            --max-maf 0.95 \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --maxDP 200 \
            --minDP 4 \
            --minQ 50 \
	    --thin 261 \
            --recode \
            --recode-INFO-all | \
      gawk -f $ADD_FLAG_DIR/add_flag.awk > flagged.vcf
   fi
   # It would be faster to pipe the output of add_flag.awk script to the next
   # script. But I want to keep the flagged.vcf file, in case we need a different
   # filtering strategy. The name of the output file encodes the minimum number
   # of individuals requested in each species. E.g., r10e10c5.vcf means 10 roumanicus,
   # 10 europaeus, and 5 concolor.
   gawk -f filtervcf.awk popmap.txt flagged.vcf > r10e10c4.vcf
fi
