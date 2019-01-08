#!/bin/bash
#
#				2018-12-15
#				----------
#
# In previous folders, it became clear that introgression lingers in some genomic
# regions longer than in others. However, it would be impossible to determine what
# alleles from one ancestry are more likely to keep segregating in the other ancesty's
# background, using only one admixed individual. Barbora had the very good idea of
# estimating genomic clines. Genomic clines describe how probable it is for a specific
# SNP to have an alternative ancestry, as a function of the overall genomic proportion
# of introgressed variation. Z. Gompert and C.A. Buerkle (2012, Mol. Ecol. Res. doi:
# 10.1111/1755-0998.12009.x), devised a Bayesian method to estimate genomic clines,
# which uses two ancestral populations and one or more admixed populations.
#
# First, we need to define a representative sample of E. europaeus with the least possible
# levels of introgression, and a similarly 'pure' sample of E. roumanicus. Other specimens,
# likely those in the hybrid zone and including the known admixed ones, would represent
# different levels of admixture. The program estimates simultaneously the admixture index
# of each individual and the parameters that describe the clines in each SNP. Note that
# there is a large number of parameters to be estimated (at least 2 per SNP). This method
# gains power from larger number of individuals, rather than from longer lists of SNPs.
# actually, for computational reasons, it may be necessary to select a subset of the
# SNPs. In that case, it makes sense to focus on the SNPs that better differentiate the
# two ancestries, namely, those that are fixed for different alleles in the two species.
#
#
# Definition of populations
# -------------------------
#
# In principle, geography could inform us of what individuals should be more pure breeds
# and what others have a chance of being partially admixed (hybrid zone). However, it is
# clear from the Admixture analysis that recent admixture in the hybrid zone is rare.
# And the abba/baba tests suggested that some introgression between E. roumanicus and E.
# europaeus must be older than the present interglacial period, since western E. europaeus
# also shows signs of introgression. Thus, it is unlikely that geography alone can tell
# us what individuals carry any introgressed variation.
#
# An alternative approach is to perform a PCA analysis with only E. roumanicus and E. europaeus
# species, so that the main axis of variation should represent the genetic positions of
# the samples between the two pure lineages. However, this approach can be confused by
# the presence of additional genetic contributions. For example, individuals from islands
# and individuals with a higher proportion of missing data may not be correctly placed within
# that axis of genetic variation.
#
# I think that a combined approach would be better. There is a trade off between the number of
# loci that can be used (to place individuals in the correct group) and the number of individuals
# that can be used in the estimation of the genomic clines. Note that genomic clines are
# estimated for every SNP. Thus, power is not increased by using more SNPs, but by including
# more individuals. Plus, the Bayesian analysis may be prohibitively long with too many SNPs.

GZVCF=/data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz
POPDATA=../../data/populations.txt

if [ ! -e popmap.txt ]; then
   gawk 'BEGIN{
      POPNAME[1] = "roumanicus"
      POPNAME[2] = "europaeus"
      POPNAME[3] = "concolor"
      POPNAME[4] = "hybrid"
      POPNAME[5] = "Hemiechinus"
      POPNAME[6] = "Atelerix"
   }(/^[^#]/){
      print $1 "\t" POPNAME[$2]
   }' $POPDATA | sort -k 2,2 > popmap.txt
fi

if [ ! -e excluded.txt ] || [ $(cat excluded.txt | wc -l) -gt 11 ]; then
   echo Er65_IS25      > excluded.txt    # Hemiechinus
   echo Er73_SNG1     >> excluded.txt    # Atelerix
   echo Er59_FR1      >> excluded.txt    # wrong
   echo Er68_PRT1B    >> excluded.txt    # wrong
   echo Er38_LB1      >> excluded.txt    # concolor
   echo Er49_GR38     >> excluded.txt    # concolor
   echo Er53_ASR7     >> excluded.txt    # concolor
   echo Er64_IS1      >> excluded.txt    # concolor
   echo Er75_TRC2A    >> excluded.txt    # concolor
   echo Er27_SK32     >> excluded.txt    # admixed concolor
   echo Er26_JUG4     >> excluded.txt    # admixed concolor
fi

MAX_MISSING=0.05
if [ ! -e summary_missingness.txt ]; then
   echo -e "# Most_missing: Largest proportion of missing data per individual before the exclusion of the following specimen. The first 1.0 is not true." > summary_missingness.txt
   echo -e "# Excluding: Sample excluded in the current round." >> summary_missingness.txt
   echo -e "# Roumanicus: Number of E. roumanicus left after exclusion." >> summary_missingness.txt
   echo -e "# Europaeus: Number of E. europaeus left after exclusion." >> summary_missingness.txt
   echo -e "# Others: Number of other specimens of other species or populations left after exclusion. Usually, the one hybrid." >> summary_missingness.txt
   echo -e "# Num_sites: Number of sites with no more than two missing genotypes, after excluding the sample shown on column 2." >> summary_missingness.txt
   echo -e "Most_missing\tExcluding\tRoumanicus\tEuropaeus\tOthers\tNum_sites" >> summary_missingness.txt

   if [ ! -e EuropaeusRoumanicus.recode.vcf ]; then
      vcftools --gzvcf $GZVCF \
               --out EuropaeusRoumanicus \
               --remove excluded.txt \
               --maf 0.01 \
               --max-maf 0.95 \
               --remove-indels \
               --min-alleles 2 \
               --max-alleles 2 \
               --maxDP 200 \
               --minDP 4 \
               --minQ 50 \
               --thin 261 \
               --max-missing-count 2 \
               --recode \
               --recode-INFO-all
   fi
   gawk '(/^#CHROM/){
      for (i = 10; i <= NF; i++) {
         SAMPLE[i] = $i
      }
   }(/^[^#]/){
      for (i = 10; i <= NF; i++) {
         if ($i ~ /^\.$|\.\/\./) {
            MISSING[SAMPLE[i]] += 1
         }
      }
      TOTAL += 1
   }END{
      for (i in SAMPLE) {
         printf("%s\t%.4f\n", SAMPLE[i], MISSING[SAMPLE[i]] / TOTAL)
      }
   }' EuropaeusRoumanicus.recode.vcf | sort -nrk 2,2 > missingness.txt

   CURRENT_MISSING=$(head -n 1 missingness.txt | cut -f 2)

   # Note the use of bc to compare floats. Using double parentheses is equivalent to and shorter than:
   # while [ $(echo "$CURRENT_MISSING>$MAX_MISSING" | bc -l) -eq 1 ]; do
   while (( $(echo "$CURRENT_MISSING>$MAX_MISSING" | bc -l) )); do
      EXCLUDE=$(head -n 1 missingness.txt | cut -f 1)
      mv missingness.txt missingness_$CURRENT_MISSING.txt
      echo $EXCLUDE >> excluded.txt
      if [ ! -e EuropaeusRoumanicus_$CURRENT_MISSING.recode.vcf ]; then
         vcftools --gzvcf $GZVCF \
                  --out EuropaeusRoumanicus_$CURRENT_MISSING \
                  --remove excluded.txt \
                  --maf 0.01 \
                  --max-maf 0.95 \
                  --remove-indels \
                  --min-alleles 2 \
                  --max-alleles 2 \
                  --maxDP 200 \
                  --minDP 4 \
                  --minQ 50 \
                  --thin 261 \
                  --max-missing-count 2 \
                  --recode \
                  --recode-INFO-all
      fi

      gawk '(/^#CHROM/){
         for (i = 10; i <= NF; i++) {
            SAMPLE[i] = $i
         }
      }(/^[^#]/){
         for (i = 10; i <= NF; i++) {
            if ($i ~ /^\.$|\.\/\./) {
               MISSING[SAMPLE[i]] += 1
            }
         }
         TOTAL += 1
      }END{
         for (i in SAMPLE) {
            printf("%s\t%.4f\n", SAMPLE[i], MISSING[SAMPLE[i]] / TOTAL)
         }
      }' EuropaeusRoumanicus_$CURRENT_MISSING.recode.vcf | sort -nrk 2,2 > missingness.txt

      ROUMANICUS=$(cut -f 1 missingness.txt | grep -f - popmap.txt | grep roumanicus | wc -l)
      EUROPAEUS=$(cut -f 1 missingness.txt | grep -f - popmap.txt | grep europaeus  | wc -l)
      OTHERS=$(cut -f 1 missingness.txt | grep -f - popmap.txt | grep -vP "roumanicus|europaeus" | wc -l)
      SITES=$(grep -v "^#" EuropaeusRoumanicus_$CURRENT_MISSING.recode.vcf | wc -l)
      echo -e "$CURRENT_MISSING\t$EXCLUDE\t$ROUMANICUS\t$EUROPAEUS\t$OTHERS\t$SITES" >> summary_missingness.txt

      if [ ! -e scores_$CURRENT_MISSING.txt ]; then
         if [ ! -e EuropaeusRoumanicus_$CURRENT_MISSING.012 ]; then
            vcftools --vcf EuropaeusRoumanicus_$CURRENT_MISSING.recode.vcf --012 --out EuropaeusRoumanicus_$CURRENT_MISSING
         fi
         sed -i 's/-1/NA/g' EuropaeusRoumanicus_$CURRENT_MISSING.012
         R --save < run_Adegenet.R --args EuropaeusRoumanicus_$CURRENT_MISSING.012 $SITES EuropaeusRoumanicus_$CURRENT_MISSING.012.indv EuropaeusRoumanicus_$CURRENT_MISSING.012.pos scores_$CURRENT_MISSING.txt
         LC_ALL=C sort -gk 2,2 scores_$CURRENT_MISSING.txt > z1
         mv z1 scores_$CURRENT_MISSING.txt
      fi
      CURRENT_MISSING=$(head -n 1 missingness.txt | cut -f 2)
   done
fi

if [ -e missingness.txt ]; then rm missingness.txt; fi
if [ -e excluded.txt ];    then rm excluded.txt;    fi

