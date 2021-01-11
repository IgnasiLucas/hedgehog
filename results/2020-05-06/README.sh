#!/bin/bash
#
#				2020-05-06
#				==========
#
# The goal here is to reproduce Kristyna's analysis from 2020-04-13. It is an
# application of Bayescan to identify candidate loci under natural selection
# among populations of hedgehogs.
#
# ============================================================================
# DATA
# ============================================================================
#
# I copy the following files in a ../../data/2020-05-06 folder, to make sure
# they are available for subsequent runs: popmap_all.txt and all_merged_2020.vcf.

DATADIR="../../data/2020-05-06"
KRISTYNA="/data/kristyna/hedgehog/results_2020/2020-04-13"

# ============================================================================
# FILTERING VCF
# ============================================================================
#
# Only E. europaeus, E. roumanicus and admixed individuals should be present in
# the vcf file. I will use an "include" file, instead of a "remove" one.

if [ ! -e era.recode.vcf ]; then
   if [ ! -e include.txt ]; then
      gawk '(($1 ~ /^Er/) && ($2 ~ /[127]/)){print $1}' $DATADIR/popmap_all.txt > include.txt
   fi

   if [ ! -e first.recode.vcf ]; then
      vcftools --vcf $DATADIR/all_merged_2020.vcf \
               --keep include.txt \
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
               --out first
   fi

   if [ ! -e first.imiss ]; then
      vcftools --vcf first.recode.vcf --missing-indv --out first
   fi

   if [ ! -e exclude_incomplete.txt ]; then
      gawk '((NR > 1) && ($5 + 0 >= 0.6)){print $1}' first.imiss > exclude_incomplete.txt
   fi

   vcftools --vcf $DATADIR/all_merged_2020.vcf \
            --keep include.txt \
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
            --out era

   vcftools --vcf era.recode.vcf --missing-indv --out era
fi

# =============================================================================
# Preparing input files for bayescan
# =============================================================================
#
# BayeScan2.1 is now installed in /usr/local/bin, available to all users by typing
# "bayescan". It requires SNP data in a fancy format. I paste below a sample of it:
#
#     [loci]=1000
#
#     [populations]=10
#
#     [pop]=1
#       1  40  2    0  40
#       2  40  2   34  6
#       3  40  2   19  21
#       4  40  2   31  9
#       5  40  2   14  26
#       6  40  2   40  0
#       7  40  2   29  11
#       8  40  2   20  20
#       9  40  2    9  31
#      10  40  2   32  8
#      11  40  2   26  14
#      12  40  2   17  23
#      13  40  2   27  13
#
# First two lines show number of loci and number of populations. Then, for each
# population, one locus per line, with: locus number, total number of sampled
# genes (twice the number of individuals with data if diploids), number of alleles
# (2: reference and alternative), absolute frequency of each allele.
#
# I will try to generate the files here, with awk. For simplicity, I will read the
# vcf file three times, once per population, and then join the pairs in two different
# input files: one with E. roumanicus and admixed individuals, and one with E.
# europaeus and the admixed individuals.

for pop in europaeus roumanicus admixed; do
   if [ ! -e $pop.bayescan.txt ]; then
      gawk -v POP=$pop 'BEGIN{
         POPNUM["europaeus"] = 1
         POPNUM["roumanicus"] = 2
         POPNUM["admixed"] = 7
         LOCUS = 1
      }((FILENAME ~ /popmap_all/) && (/^Er/)){
         if ($2 == POPNUM[POP]) MEMBER[$1] = 1
      }((FILENAME ~ /vcf$/) && (/^#CHROM/)){
         for (i = 10; i <= NF; i++) {
            if ($i in MEMBER) INCLUDE[i] = 1
         }
      }((FILENAME ~ /vcf$/) && (/^[^#]/)){
         REF = 0
         ALT = 0
         for (i in INCLUDE) {
            if ($i ~ /^0\/0/) REF += 2
            if ($i ~ /^0\/1/) {REF += 1; ALT++}
            if ($i ~ /^1\/1/) ALT += 2
         }
         print LOCUS "\t" REF + ALT "\t2\t" REF "\t" ALT
         LOCUS++
      }' $DATADIR/popmap_all.txt era.recode.vcf > $pop.bayescan.txt
   fi
done

# The code above seems to produce valide output. Comparing it to Kristyna's output
# from a program called PGDspider, I note that PGDspider orders alleles alphabetically,
# while I report always the reference allele first, and then the alternative. Another
# difference is that PGDspider reports monomorphic loci as having only one allele. To
# do that, I should check that all populations are monomorphic. I think BayeScan should
# tolerate my output anyways.

if [ ! -d europaeus_admixed ];  then mkdir europaeus_admixed;  fi
if [ ! -d roumanicus_admixed ]; then mkdir roumanicus_admixed; fi
if [ ! -d roumanicus_europaeus ]; then mkdir roumanicus_europaeus; fi
if [ ! -d roumanicus_europaeus_admixed ]; then mkdir roumanicus_europaeus_admixed; fi

NUMLOCI=$(tail -n 1 europaeus.bayescan.txt | cut -f 1)

if [ ! -e roumanicus_admixed/ra_bayescan.txt ]; then
   echo -e "[loci]=$NUMLOCI\n" >  roumanicus_admixed/ra_bayescan.txt
   echo -e "[populations]=2\n" >> roumanicus_admixed/ra_bayescan.txt
   echo -e "[pop]=1"           >> roumanicus_admixed/ra_bayescan.txt
   cat roumanicus.bayescan.txt >> roumanicus_admixed/ra_bayescan.txt
   echo -e "\n[pop]=2"         >> roumanicus_admixed/ra_bayescan.txt
   cat admixed.bayescan.txt    >> roumanicus_admixed/ra_bayescan.txt
fi

if [ ! -e europaeus_admixed/ea_bayescan.txt ]; then
   echo -e "[loci]=$NUMLOCI\n" >  europaeus_admixed/ea_bayescan.txt
   echo -e "[populations]=2\n" >> europaeus_admixed/ea_bayescan.txt
   echo -e "[pop]=1"           >> europaeus_admixed/ea_bayescan.txt
   cat europaeus.bayescan.txt  >> europaeus_admixed/ea_bayescan.txt
   echo -e "\n[pop]=2"         >> europaeus_admixed/ea_bayescan.txt
   cat admixed.bayescan.txt    >> europaeus_admixed/ea_bayescan.txt
fi

if [ ! -e roumanicus_europaeus/re_bayescan.txt ]; then
   echo -e "[loci]=$NUMLOCI\n"  > roumanicus_europaeus/re_bayescan.txt
   echo -e "[populations]=2\n" >> roumanicus_europaeus/re_bayescan.txt
   echo -e "[pop]=1"           >> roumanicus_europaeus/re_bayescan.txt
   cat roumanicus.bayescan.txt >> roumanicus_europaeus/re_bayescan.txt
   echo -e "\n[pop]=2"         >> roumanicus_europaeus/re_bayescan.txt
   cat europaeus.bayescan.txt  >> roumanicus_europaeus/re_bayescan.txt
fi

if [ ! -e roumanicus_europaeus_admixed/rea_bayescan.txt ]; then
   echo -e "[loci]=$NUMLOCI\n"  > roumanicus_europaeus_admixed/rea_bayescan.txt
   echo -e "[populations]=3\n" >> roumanicus_europaeus_admixed/rea_bayescan.txt
   echo -e "[pop]=1"           >> roumanicus_europaeus_admixed/rea_bayescan.txt
   cat roumanicus.bayescan.txt >> roumanicus_europaeus_admixed/rea_bayescan.txt
   echo -e "\n[pop]=2"         >> roumanicus_europaeus_admixed/rea_bayescan.txt
   cat europaeus.bayescan.txt  >> roumanicus_europaeus_admixed/rea_bayescan.txt
   echo -e "\n[pop]=3"         >> roumanicus_europaeus_admixed/rea_bayescan.txt
   cat admixed.bayescan.txt    >> roumanicus_europaeus_admixed/rea_bayescan.txt
fi

# =====================================================================================
# Running BayeScan2.1
# =====================================================================================

if [ ! -e europaeus_admixed/ea_fst.txt ]; then
   bayescan europaeus_admixed/ea_bayescan.txt \
            -threads 8 \
            -o ea \
            -od ./europaeus_admixed \
            -n 50000 \
            -thin 50 \
            -nbp 20 \
            -pilot 5000 \
            -burn 50000 \
            -all_trace \
            -pr_odds 10
fi

if [ ! -e roumanicus_admixed/ra_fst.txt ]; then
   bayescan roumanicus_admixed/ra_bayescan.txt \
            -threads 8 \
            -o ra \
            -od ./roumanicus_admixed \
            -n 50000 \
            -thin 50 \
            -nbp 20 \
            -pilot 5000 \
            -burn 50000 \
            -all_trace \
            -pr_odds 10
fi

#if [ ! -e roumanicus_europaeus/re_fst.txt ]; then
#   bayescan roumanicus_europaeus/re_bayescan.txt \
#            -threads 8 \
#            -o ra \
#            -od ./roumanicus_europaeus \
#            -n 50000 \
#            -thin 50 \
#            -nbp 20 \
#            -pilot 5000 \
#            -burn 50000 \
#            -all_trace \
#            -pr_odds 10
#fi
#
#if [ ! -e roumanicus_europaeus_admixed/rea_fst.txt ]; then
#   bayescan roumanicus_europaeus_admixed/rea_bayescan.txt \
#            -threads 8 \
#            -o ra \
#            -od ./roumanicus_europaeus_admixed \
#            -n 50000 \
#            -thin 50 \
#            -nbp 20 \
#            -pilot 5000 \
#            -burn 50000 \
#            -all_trace \
#            -pr_odds 10
#fi
#
# =============================================================================
# COMPARISON WITH KRISTYNA'S RESULTS
# =============================================================================
#
# Before getting to understand BayeScan's output any deeper, I want to accomplish
# the original purpose of this script, which was merely to check Kristyna's work
# from 2020-04-13. I will compare the "*_fst.txt" outputs in R and I will make some
# plots to make sense of the results.

if [ ! -e loci.txt ]; then
   gawk '(/^[^#]/){print $1 "\t" $2}' era.recode.vcf > loci.txt
fi

if [ ! -e bayescan.html ]; then
   R -q --no-save -e "renderIt <- function(kea, kra, kera, iea, ira, iera) {
                         rmarkdown::render('bayescan.Rmd',
                         params = list(KEA=kea, KRA=kra, KERA=kera, IEA=iea, IRA=ira, IERA=iera),
                         output_file = 'bayescan.html')
                      }" \
     -e "renderIt(kea  = '$KRISTYNA/europaeus_admixed_bayescan10_fst.txt',
                  kra  = '$KRISTYNA/roumanicus_admixed_bayescan10_fst.txt',
                  kera = '$KRISTYNA/era_bayescan10_fst.txt',
                  iea  = 'europaeus_admixed/ea_fst.txt',
                  ira  = 'roumanicus_admixed/ra_fst.txt',
                  iera = 'roumanicus_europaeus_admixed/rea_fst.txt')"
fi

# =============================================================================
# CONCLUSIONS
# =============================================================================

