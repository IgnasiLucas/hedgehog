#!/bin/bash
#
#				2018-07-17
#				==========
#
# The Admixture analysis used a set of sites selected from the original vcf
# file with the following filters:
#
#    vcftools --gzvcf $SOURCE/merged.vcf.gz \
#            --keep popmap3 \
#            --recode \
#            --maf 0.01 \
#            --max-maf 0.95 \
#            --remove-indels \
#            --min-alleles 2 \
#            --max-alleles 2 \
#            --maxDP 200 \
#            --minDP 6 \
#            --minGQ 20 \
#            --max-missing 0.8 \
#            --recode-INFO-all \
#            --out $DATADIR/erinaceus
#
# (Kristyna, 05-06-2018). This are relatively stringent filters meant to
# obtain a reasonably complete matrix of genotypes. It produced 67459 sites.
# These are the sites that I used later to identify the 2033 SNPs with
# europaeus ancestry in the hybrid's genome (2018-06-30).
#
# On the other side, the abba/baba test did not require a complete matrix
# of genotypes, and allowed us to use less stringent filters. In addition,
# I used a bed file with genomic regions of well-behaved total coverage to
# exclude highly repetitive or unreliable parts of the reference genome.
# The exact filters used to prepare the abba/baba test were this:
#
#    vcftools --vcf flagged.vcf \
#            --out thinned \
#            --bed goodregions.bed \
#            --minQ 50.0 \
#            --thin 261 \
#            --recode \
#            --recode-INFO NS \
#            --recode-INFO AN \
#            --recode-INFO AF \
#            --recode-INFO ABP \
#            --recode-INFO BPF
#
# It produced 668416 sites, which I hoped would include the 2033 sites
# mentioned before. All these sites were further filtered by our script
# freq2.awk to select sites appropriate for abba/baba test, with a minimum
# number of genotypes from each population, and with a presumed ancestral
# allele fixed in the outgroup. The exact number of sites actually used
# depended on the definition of populations P1, P2 and P3: between 180279
# and 369451. Kristyna (17-07-2018) noticed that the number of sites identified
# as having europaeus ancestry in the hybrid that were included in the abba/baba
# test was very low: 83.
#
# The first question is what filters removed them.

ABBABABA_VCF=/data/joiglu/hedgehog/results/2018-03-27b/thinned.recode.vcf
ADMIXTURE_VCF=/data/kristyna/hedgehog/results_2018/05-06-2018/erinaceus.recode.vcf
GOODREGIONS=../2018-03-27/pooled_filtered.bed
ORIGINAL_VCF=/data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz
BIN=../../bin

if [ ! -e 2033.txt ]; then
   gawk '{print $1 "\t" $2 "\t"}' ../2018-06-30/Er37_SK27.europaeus_SNPs.txt > 2033.txt
fi

if [ ! -e summary.txt ]; then
   if [ ! -e in_vcf.txt ]; then
      grep -Ff 2033.txt $ABBABABA_VCF | gawk '{print $1 "\t" $2 "\t" }' > in_vcf.txt
   fi
   IN_VCF=`cat in_vcf.txt | wc -l`
   if [ ! -e missing.txt ]; then
      grep -Fvf in_vcf.txt 2033.txt > missing.txt
   fi
   if [ ! -e low_Q.txt ]; then
      grep -Ff missing.txt $ADMIXTURE_VCF | gawk '($6 < 50){print $1 "\t" $2 "\t"}' > low_Q.txt
   fi
   LOW_Q=`cat low_Q.txt | wc -l`
   if [ ! -e missing_highQ.txt ]; then
      grep -Fvf low_Q.txt missing.txt > missing_highQ.txt
   fi
   if [ ! -e missing_highQ.bed ]; then
      gawk '{print $1 "\t" $2 "\t" $2}' missing_highQ.txt > missing_highQ.bed
   fi
   if [ ! -e bad_regions.txt ]; then
      bedtools intersect -wa -v -a missing_highQ.bed -b $GOODREGIONS | gawk '{print $1 "\t" $2 "\t"}' > bad_regions.txt
   fi
   BAD_REGIONS=`cat bad_regions.txt | wc -l`
   if [ ! -e missing_highQ_goodregions.txt ]; then
      grep -Fvf bad_regions.txt missing_highQ.txt > missing_highQ_goodregions.txt
   fi
   if [ ! -e missing_highQ_goodregions_thinnedout.txt ]; then
      if [ ! -e thick.recode.vcf ]; then
         if [ ! -e goodregions.bed ]; then
            echo "track name=\"bamcoverage\" description=\"Fragments covered between 200 and 1000 times among all 50 samples\"" > goodregions.bed
            cat $GOODREGIONS >> goodregions.bed
         fi
         if [ ! -e flagged.vcf ]; then
            gunzip -c $ORIGINAL_VCF > z1.vcf
            gawk -f add_flag.awk z1.vcf > flagged.vcf
            rm z1.vcf
         fi
         vcftools --vcf flagged.vcf \
               --out thick \
               --bed goodregions.bed \
               --minQ 50.0 \
               --recode \
               --recode-INFO NS \
               --recode-INFO AN \
               --recode-INFO AF \
               --recode-INFO ABP \
               --recode-INFO BPF
      fi
      grep -Ff missing_highQ_goodregions.txt thick.recode.vcf | gawk '{print $1 "\t" $2 "\t"}' > missing_highQ_goodregions_thinnedout.txt
   fi
   THINNED=`cat missing_highQ_goodregions_thinnedout.txt | wc -l`

   printf "%10d sites removed because of low quality.\n" $LOW_Q > summary.txt
   printf "%10d sites removed because in regions with abnormal coverage.\n" $BAD_REGIONS >> summary.txt
   printf "%10d sites removed while thinning the vcf out to reduce linkage disequilibrium issues.\n" $THINNED >> summary.txt
   printf "%10d sites were left in the vcf for the abba/baba test.\n" $IN_VCF >> summary.txt
   printf "  ---------\n" >> summary.txt
   printf "%10d total, which should be equal to 2033 sites with europaeus ancestry in hybrid.\n" $(( $LOW_Q + $BAD_REGIONS + $THINNED + $IN_VCF )) >> summary.txt
fi

# CONCLUSIONS
# ===========
#
# These are the results.
#
#        64 sites removed because of low quality.
#       297 sites removed because in regions with abnormal coverage.
#      1512 sites removed while thinning the vcf out to reduce linkage disequilibrium issues.
#       160 sites were left in the vcf for the abba/baba test.
#  ---------
#      2033 total, which should be equal to 2033 sites with europaeus ancestry in hybrid.
#
#
# After talking with Kristyna, we decided to re-do the filtering of the original vcf file
# and repeat all downstream analyses from the newly filtered vcf. She applied the following
# filters in 24-07-2018:
#
#   vcftools --gzvcf $SOURCE/merged.vcf.gz \
#            --keep popmap \
#            --recode \
#            --maf 0.01 \
#            --max-maf 0.95 \
#            --remove-indels \
#            --min-alleles 2 \
#            --max-alleles 2 \
#            --maxDP 200 \
#            --minDP 4 \
#            --minQ 50 \
#	    --thin 261 \
#            --max-missing 0.75 \
#            --recode-INFO-all \
#            --out $DATADIR/all
#

