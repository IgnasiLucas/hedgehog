#!/bin/bash
#
#				2018-03-27
#				----------
#
# We received more sequence data, and Kristýna has already aligned it and
# merged it with the previous data. Here, I want to compare the completeness
# of the vcf files before and after the last sequencing run.
#
# The last vcf build before the last sequencing run is ../2017-05-25/hhog.vcf.
# It was obtained with freebayes. It includes a binary flag indicating what
# samples each allele is present in. It is not filtered by quality.

BEFORE=../2017-05-25/hhog.vcf

# The vcf obtained by Kristýna after the addition of the last batch of sequences
# is in /data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz. It was
# also obtained with freebayes, with similar options.

AFTER=/data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz

# For the comparison to be accurate and fast, I could use only a subset of sites.
# For example, the sites covered in the only Hemichinus sample, which are the
# sites among which the abba/baba test is feasible.

if [ ! -e popmap.txt ]; then
   cp ../2017-05-25/populations.txt ./popmap.txt
fi

# The order of samples in the vcf files is not the same. I will select the sites
# where Hemichinus was covered in the vcf file before the last sequencing event,
# and the quality of the variant was at least 50.

HEMICHINUS_POSITION_BEFORE=$( head -n 6000 $BEFORE | gawk '(/^#CHROM/){for (i = 0; i <= NF - 10; i++) {if ($(i + 10) == "Er65_IS25") print i}}' )
echo "Hemichinus is sample number $HEMICHINUS_POSITION_BEFORE"
HEMICHINUS_FLAG=$( echo "2 ^ $HEMICHINUS_POSITION_BEFORE" | bc )
echo "The corresponding binary number is $HEMICHINUS_FLAG"

if [ ! -e comparison.txt ]; then
   if [ ! -e hemichinus_sites.txt ]; then
      gawk -f filtervcf.awk -v REQUIRED=$HEMICHINUS_FLAG $BEFORE | gawk '((/^[^#]/) && ($6 >= 50.0)){print $1 "\t" $2}' > hemichinus_sites.txt
   fi

   if [ ! -e before.recode.vcf ] && [ ! -e common_before.recode.vcf ]; then
      vcftools --vcf $BEFORE \
               --out before \
               --positions hemichinus_sites.txt \
               --recode 1> before.log 2> before.err &
   fi

   if [ ! -e after.recode.vcf ] && [ ! -e common_after.recode.vcf ]; then
      vcftools --gzvcf $AFTER \
               --out after \
               --positions hemichinus_sites.txt \
               --recode 1> after.log 2> after.err &
   fi

   wait

   # In order to use exactly the same sites, I will reduce the previous files
   # further. First, I need to exclude chromosomes present in only one file.
   # I manually determined them and use them below.

   if [ ! -e beforeafter.diff.sites_in_files ]; then
      vcftools --vcf before.recode.vcf \
               --out beforeafter \
               --not-chr NW_006805355.1 \
               --not-chr NW_006805505.1 \
               --not-chr NW_006805518.1 \
               --not-chr NW_006805635.1 \
               --not-chr NW_006805650.1 \
               --not-chr NW_006805660.1 \
               --not-chr NW_006805683.1 \
               --not-chr NW_006805778.1 \
               --not-chr NW_006805784.1 \
               --not-chr NW_006805791.1 \
               --not-chr NW_006805810.1 \
               --not-chr NW_006805815.1 \
               --not-chr NW_006805832.1 \
               --not-chr NW_006805936.1 \
               --not-chr NW_006806027.1 \
               --not-chr NW_006806028.1 \
               --not-chr NW_006806299.1 \
               --not-chr NW_006806372.1 \
               --not-chr NW_006806396.1 \
               --not-chr NW_006806515.1 \
               --not-chr NW_006806796.1 \
               --not-chr NW_006807000.1 \
               --not-chr NW_006807010.1 \
               --not-chr NW_006807104.1 \
               --not-chr NW_006807718.1 \
               --not-chr NW_006807727.1 \
               --not-chr NW_006807972.1 \
               --not-chr NW_006807985.1 \
               --not-chr NW_006808106.1 \
               --not-chr NW_006808186.1 \
               --not-chr NW_006808370.1 \
               --not-chr NW_006809027.1 \
               --not-chr NW_006809075.1 \
               --not-chr NW_006809618.1 \
               --diff after.recode.vcf \
               --diff-site
   fi

   if [ ! -e common_sites.txt ]; then
      grep B beforeafter.diff.sites_in_files | cut -f 1,2 > common_sites.txt
   fi

   if [ ! -e common_before.recode.vcf ]; then
      vcftools --vcf before.recode.vcf \
               --out common_before \
               --positions common_sites.txt \
               --recode 1> z1.log 2> z2.err &
      rm before.recode.vcf
   fi

   if [ ! -e common_after.recode.vcf ]; then
      vcftools --vcf after.recode.vcf \
               --out common_after \
               --positions common_sites.txt \
               --recode 1> z3.log 2> z4.err &
      rm after.recode.vcf
   fi

   wait

   # Now I want to create a file with at least six columns: chromosome,
   # site, variant quality before addition of data, variant quality
   # afterwards, number of samples with data before, and after data
   # addition.

   gawk '(/^[^#]/){
      WITHDATA = 0
      for (i = 10; i <= NF; i++) {
         if ($i ~ /[0-9]/) WITHDATA++
      }
      print $1 "\t" $2 "\t" $6 "\t" WITHDATA
   }' common_before.recode.vcf > z1 &
   gawk '(/^[^#]/){
      WITHDATA = 0
      for (i = 10; i <= NF; i++) {
         if ($i ~ /[0-9]/) WITHDATA++
      }
      print $6 "\t" WITHDATA
   }' common_after.recode.vcf > z2 &
   wait
   paste z1 z2 > comparison.txt
   rm z1 z2
fi
rm 

# Visual inspection of the comparison file, with R, shows that there
# is substantial improvement. The average number of individuals with
# data per site changed from 32.6 to 40.5. The median changed from
# 37 to 44. The qualities also improved.
#
# However, some sites have fewer samples with data and/or lower qualities
# after the addition of sequencing data. A list of possible reasons for that:
#
#   * The filter settings in freebayes were different to produce the two original
#     vcf files. Before the addition of sequence data, the minimum mapping quality
#     was 30 and the minimum base quality, 20. After the addition of sequence
#     data, both thresholds were set to 15. Thus, part of the increase in the
#     number of samples with data may be due to the less stringent use of the
#     available reads. That could also explain some lower qualities. But it does
#     not help explain why in some sites there are fewer samples with data.
#
#   * The mapping of reads may have been different. I need to check that.
#

