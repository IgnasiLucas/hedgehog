#!/bin/bash
#
#				2020-06-26
#				==========
#
# One difficulty with bayescan is to select the individuals in each population.
# The popmap files are not very stable, due to occasional changes in the adscription
# of individuals to populations. There should be a script somewhere that takes the
# output from structure and revises the original popmap to fix all problems. But that
# is beyond the purpose of this script. I just need to check that bayescan works as
# expected, and I will use for that purpose any popmap file that I can trust in Kristýna's
# folders.
#
# The file kristyna/results_2020/2020-06-16/popmap_all_correct seems to be reliable and
# extensively used. It includes central European and Russian samples, and has some samples
# labeled as wrong. In addition, individuals Er61_GR87 and Er62_GR95 are supposed to have
# too much missing data, but I will take care of that later. I copy that file to the data
# directory.
#
# I note that the labeling of populations changed between the two main data sets: when
# we had only Central European samples, E. roumanicus was 1 and E. europaeus, 2. Now,
# they are reversed. In any case, it is clear from visual inspection of the popmap files
# that the different labelings were never combined. Clearly, "popmap_all_correct.txt"
# must be correct.

POPMAP='../../data/2020-06-26/popmap_all_correct.txt'
VCF='/data/kristyna/hedgehog/results_2020/2020-02-28/all_merged_2020.vcf'

# The popmap file does not include the bioregion of origin (central Europe or Russia).
# Kristýna assumes samples numbered from 76 on are Russian. Popmap files in
# kristyna/results_2020/2020-05-02 confirms it.

if [ ! -e CE_europaeus.txt ]; then
   gawk '(/^Er/){
      match($1, /Er([0-9]+)_/, A)
      if ((A[1] + 0 < 76) && ($2 == "1")) print $1
   }' $POPMAP > CE_europaeus.txt
fi

if [ ! -e CE_roumanicus.txt ]; then
   gawk '(/^Er/){
      match($1, /Er([0-9]+)_/, A)
      if ((A[1] + 0 < 76) && ($2 == "2")) print $1
   }' $POPMAP > CE_roumanicus.txt
fi

if [ ! -e RB_europaeus.txt ]; then
   gawk '(/^Er/){
      match($1, /Er([0-9]+)_/, A)
      if ((A[1] + 0 >= 76) && ($2 == "1")) print $1
   }' $POPMAP > RB_europaeus.txt
fi

if [ ! -e RB_roumanicus.txt ]; then
   gawk '(/^Er/){
      match($1, /Er([0-9]+)_/, A)
      if ((A[1] + 0 >= 76) && ($2 == "2")) print $1
   }' $POPMAP > RB_roumanicus.txt
fi

#if [ ! -e RB_admixed.txt ]; then
#   gawk '(/^Er/){
#      match($1, /Er([0-9]+)_/, A)
#      if ((A[1] + 0 >= 76) && ($2 == "7")) print $1
#   }' $POPMAP > RB_admixed.txt
#fi
#
#if [ ! -e admixed.txt ]; then
#   gawk '(/^Er/){
#       if ($2 == "7") print $1
#   }' $POPMAP > admixed.txt
#fi

# =================================================================
# FILTERING VCF
# =================================================================
#
# I run vcftools twice: first to select the valid individuals and to
# see for how many loci individuals miss genotypes, and then to remove
# the individuals with too many missing genotypes. I do it separately
# for roumanicus and europaeus.

if [ ! -e europaeus.recode.vcf ]; then
   if [ ! -e tmp.recode.vcf ]; then
      vcftools --vcf $VCF \
               --keep CE_europaeus.txt \
               --keep RB_europaeus.txt \
               --recode \
               --maf 0.0125 \
               --max-maf 0.9875 \
               --remove-indels \
               --min-alleles 2 \
               --max-alleles 2 \
               --maxDP 200 \
               --minDP 6 \
               --minQ 50 \
               --thin 261 \
               --max-missing 0.80 \
               --recode-INFO-all \
               --out tmp
   fi

   if [ ! -e tmp.imiss ]; then
      vcftools --vcf tmp.recode.vcf --missing-indv --out tmp
   fi

   if [ ! -e exclude_incomplete.txt ]; then
      gawk '((NR > 1) && ($5 + 0 >= 0.7)){print $1}' tmp.imiss > exclude_incomplete.txt
   fi

   vcftools --vcf $VCF \
            --keep CE_europaeus.txt \
            --keep RB_europaeus.txt \
            --remove exclude_incomplete.txt \
            --recode \
            --maf 0.0125 \
            --max-maf 0.9875 \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --maxDP 200 \
            --minDP 6 \
            --minQ 50 \
            --thin 261 \
            --max-missing 0.90 \
            --recode-INFO-all \
            --out europaeus

   rm tmp.imiss exclude_incomplete.txt tmp.recode.vcf
fi

if [ ! -e roumanicus.recode.vcf ]; then
   if [ ! -e tmp.recode.vcf ]; then
      vcftools --vcf $VCF \
               --keep CE_roumanicus.txt \
               --keep RB_roumanicus.txt \
               --recode \
               --maf 0.0125 \
               --max-maf 0.9875 \
               --remove-indels \
               --min-alleles 2 \
               --max-alleles 2 \
               --maxDP 200 \
               --minDP 6 \
               --minQ 50 \
               --thin 261 \
               --max-missing 0.80 \
               --recode-INFO-all \
               --out tmp
   fi

   if [ ! -e tmp.imiss ]; then
      vcftools --vcf tmp.recode.vcf --missing-indv --out tmp
   fi

   if [ ! -e exclude_incomplete.txt ]; then
      gawk '((NR > 1) && ($5 + 0 >= 0.7)){print $1}' tmp.imiss > exclude_incomplete.txt
   fi

   vcftools --vcf $VCF \
            --keep CE_roumanicus.txt \
            --keep RB_roumanicus.txt \
            --remove exclude_incomplete.txt \
            --recode \
            --maf 0.0125 \
            --max-maf 0.9875 \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --maxDP 200 \
            --minDP 6 \
            --minQ 50 \
            --thin 261 \
            --max-missing 0.90 \
            --recode-INFO-all \
            --out roumanicus

   rm tmp.imiss exclude_incomplete.txt tmp.recode.vcf
fi

# ==================================================================================
# BAYESCAN INPUT
# ==================================================================================

if [ ! -d europaeus ]; then mkdir europaeus; fi
if [ ! -d roumanicus ]; then mkdir roumanicus; fi

for pop in europaeus roumanicus; do
   for region in CE RB; do
      if [ ! -e $pop/$region.bayescan.txt ]; then
         gawk '(FILENAME ~ /.txt/){
            MEMBER[$1] = 1
         }((FILENAME ~ /vcf/) && (/^#CHROM/)){
            LOCUS = 1
            for (i = 10; i <= NF; i++) {
               if ($i in MEMBER) INCLUDE[i] = 1
            }
         }((FILENAME ~ /vcf/) && (/^[^#]/)){
            REF = 0
            ALT = 0
            for (i in INCLUDE) {
               if ($i ~ /^0\/0/) REF += 2
               if ($i ~ /^0\/1/) {REF += 1; ALT++}
               if ($i ~ /^1\/1/) ALT += 2
            }
            print LOCUS "\t" REF + ALT "\t2\t" REF "\t" ALT
            LOCUS++
         }' ${region}_${pop}.txt ${pop}.recode.vcf > $pop/$region.bayescan.txt
      fi
   done
   if [ ! -e $pop/CE_RB.bayescan.txt ]; then
      NUMLOCI=$(tail -n 1 $pop/CE.bayescan.txt | cut -f 1)
      echo -e "[loci]=$NUMLOCI\n" >  $pop/CE_RB.bayescan.txt
      echo -e "[populations]=2\n" >> $pop/CE_RB.bayescan.txt
      echo -e "[pop]=1"           >> $pop/CE_RB.bayescan.txt
      cat $pop/CE.bayescan.txt    >> $pop/CE_RB.bayescan.txt
      echo -e "\n[pop]=2"         >> $pop/CE_RB.bayescan.txt
      cat $pop/RB.bayescan.txt    >> $pop/CE_RB.bayescan.txt
   fi
done

# ==================================================================================
# RUNNING BAYESCAN
# ==================================================================================

for pop in europaeus roumanicus; do
   if [ ! -e $pop/CE_RB_fst.txt ]; then
      bayescan $pop/CE_RB.bayescan.txt \
               -threads 8 \
               -o CE_RB \
               -od ./$pop \
               -n 50000 \
               -thin 50 \
               -nbp 20 \
               -pilot 5000 \
               -burn 50000 \
               -all_trace \
               -pr_odds 100
   fi
done
