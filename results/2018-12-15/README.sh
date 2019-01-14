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
GEODATA=../../data/coordinates.txt

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

NUM_SAMPLES=$(wc -l popmap.txt | cut -d " " -f 1)

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

NUM_EXCLUDED=$(wc -l excluded.txt | cut -d " " -f 1)
NUM_SAMPLES=$((NUM_SAMPLES - NUM_EXCLUDED))

# First, to evaluate the trade off between number of sites and number of individuals,
# I use the individual's portion of missing data to sequentially remove individuals
# from the vcf file, which increases the number of available sites with (almost) complete
# data. I need to allow for a number of missing genotypes per site, in order to identify
# the individuals that should be removed next, to optimize the number of sites. But I am
# actually interested in the number of sites with absolutely no missing data at all.

if [ ! -e summary_missingness.txt ]; then
   echo -e "# Excluded: Last sample removed from the dataset for containing too many missing sites." > summary_missingness.txt
   echo -e "# Num_samples: Current number of samples in the matrix of genotypes." >> summary_missingness.txt
   echo -e "# Roumanicus: Number of E. roumanicus samples left." >> summary_missingness.txt
   echo -e "# Europaeus: Number of E. europaeus samples left." >> summary_missingness.txt
   echo -e "# Others: Number of other specimens of other species or populations left. Usually, the one hybrid." >> summary_missingness.txt
   echo -e "# Most_missing: Largest proportion of missing data per individual." >> summary_missingness.txt
   echo -e "# Num_sites_2: Number of sites with two missing genotypes." >> summary_missingness.txt
   echo -e "# Num_sites_1: Number of sites with one missing genotype."  >> summary_missingness.txt
   echo -e "# Num_sites_0: Number of sites without any missing genotype." >> summary_missingness.txt
   echo -e "Excluded\tNum_samples\tRoumanicus\tEuropaeus\tOthers\tMost_missing\tNum_sites_2\tNum_sites_1\tNum_sites_0" >> summary_missingness.txt

   if [ ! -e EuropaeusRoumanicus_$NUM_SAMPLES.recode.vcf ]; then
      vcftools --gzvcf $GZVCF \
               --out EuropaeusRoumanicus_$NUM_SAMPLES \
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
               --recode
   fi

   if [ ! -e missingness_$NUM_SAMPLES.txt ]; then
      gawk '(/^#CHROM/){
         for (i = 10; i <= NF; i++) {
            SAMPLE[i] = $i
         }
      }(/^[^#]/){
         NUM_MISSING = 0
         for (i = 10; i <= NF; i++) {
            if ($i ~ /^\.$|\.\/\./) {
               MISSING[SAMPLE[i]] += 1
               NUM_MISSING += 1
            }
         }
         F[NUM_MISSING] += 1
         TOTAL += 1
      }END{
         for (i in SAMPLE) {
            printf("%s\t%.4f\n", SAMPLE[i], MISSING[SAMPLE[i]] / TOTAL)
         }
         for (i = 0; i <= 2; i++) {
            printf("#Sites missing %u genotypes:% 10u\n", i, F[i])
         }
      }' EuropaeusRoumanicus_$NUM_SAMPLES.recode.vcf | sort -nrk 2,2 > missingness_$NUM_SAMPLES.txt
   fi

   CURRENT_MISSING=$(head -n 1 missingness_$NUM_SAMPLES.txt | cut -f 2)
   ROUMANICUS=$(grep ^Er missingness_$NUM_SAMPLES.txt | cut -f 1 | grep -f - popmap.txt | grep roumanicus | wc -l)
   EUROPAEUS=$(grep ^Er missingness_$NUM_SAMPLES.txt | cut -f 1 | grep -f - popmap.txt | grep europaeus  | wc -l)
   OTHERS=$(grep ^Er missingness_$NUM_SAMPLES.txt | cut -f 1 | grep -f - popmap.txt | grep -vP "roumanicus|europaeus" | wc -l)
   SITES[0]=$(grep "Sites missing 0 genotypes:" missingness_$NUM_SAMPLES.txt | gawk '{print $5}')
   SITES[1]=$(grep "Sites missing 1 genotypes:" missingness_$NUM_SAMPLES.txt | gawk '{print $5}')
   SITES[2]=$(grep "Sites missing 2 genotypes:" missingness_$NUM_SAMPLES.txt | gawk '{print $5}')
   echo -e "  ----  \t$NUM_SAMPLES\t$ROUMANICUS\t$EUROPAEUS\t$OTHERS\t$CURRENT_MISSING\t${SITES[2]}\t${SITES[1]}\t${SITES[0]}" >> summary_missingness.txt

   if [ ! -e scores_$NUM_SAMPLES.txt ]; then
      # Vcftools seems to have a bug: the filter --max-missing-count is sometimes more conservative than i should.
      # Thus, I don't trust the number of sites in the 012 files will be the expected from the *.recode.vcf files.
      if [ -e EuropaeusRoumanicus_$NUM_SAMPLES.012.pos ]; then
         NUM_SITES=$(cat EuropaeusRoumanicus_$NUM_SAMPLES.012.pos | wc -l)
      else
         NUM_SITES=0
      fi
      i=0
      while [ $NUM_SITES -le $NUM_SAMPLES ]; do
         vcftools --vcf EuropaeusRoumanicus_$NUM_SAMPLES.recode.vcf --max-missing-count $i --012 --out EuropaeusRoumanicus_$NUM_SAMPLES
         NUM_SITES=$(cat EuropaeusRoumanicus_$NUM_SAMPLES.012.pos | wc -l)
         i=$((i+1))
      done
      sed -i 's/-1/NA/g' EuropaeusRoumanicus_$NUM_SAMPLES.012
      R --save < run_Adegenet.R --args EuropaeusRoumanicus_$NUM_SAMPLES.012 $NUM_SITES EuropaeusRoumanicus_$NUM_SAMPLES.012.indv EuropaeusRoumanicus_$NUM_SAMPLES.012.pos scores_$NUM_SAMPLES.txt
      LC_ALL=C sort -gk 2,2 scores_$NUM_SAMPLES.txt > z1
      mv z1 scores_$NUM_SAMPLES.txt
   fi


   STOP_WHEN_MISSING=0.05
   EXCLUDE=$(head -n 1 missingness_$NUM_SAMPLES.txt | cut -f 1)


   # Note the use of bc to compare floats. Using double parentheses is equivalent to and shorter than:
   # while [ $(echo "$CURRENT_MISSING>$STOP_WHEN_MISSING" | bc -l) -eq 1 ]; do
   while (( $(echo "$CURRENT_MISSING>$STOP_WHEN_MISSING" | bc -l) )); do
      echo $EXCLUDE >> excluded.txt
      NUM_SAMPLES=$((NUM_SAMPLES - 1))
      # Note that I need to exclude the next sample, together with all previously excluded ones,
      # again from the original vcf file. Otherwise, I would not increase the numbe of sites, of course.
      if [ ! -e EuropaeusRoumanicus_$NUM_SAMPLES.recode.vcf ]; then
         vcftools --gzvcf $GZVCF \
                  --out EuropaeusRoumanicus_$NUM_SAMPLES \
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
                  --recode
      fi

      if [ ! -e missingness_$NUM_SAMPLES ]; then
         gawk '(/^#CHROM/){
            for (i = 10; i <= NF; i++) {
               SAMPLE[i] = $i
            }
         }(/^[^#]/){
            NUM_MISSING = 0
            for (i = 10; i <= NF; i++) {
               if ($i ~ /^\.$|\.\/\./) {
                  MISSING[SAMPLE[i]] += 1
                  NUM_MISSING += 1
               }
            }
            TOTAL += 1
            F[NUM_MISSING] += 1
         }END{
            for (i in SAMPLE) {
               printf("%s\t%.4f\n", SAMPLE[i], MISSING[SAMPLE[i]] / TOTAL)
            }
            for (i = 0; i <= 2; i++) {
               printf("#Sites missing %u genotypes:% 10u\n", i, F[i])
            }
         }' EuropaeusRoumanicus_$NUM_SAMPLES.recode.vcf | sort -nrk 2,2 > missingness_$NUM_SAMPLES.txt
      fi

      CURRENT_MISSING=$(head -n 1 missingness_$NUM_SAMPLES.txt | cut -f 2)
      ROUMANICUS=$(grep ^Er missingness_$NUM_SAMPLES.txt | cut -f 1 | grep -f - popmap.txt | grep roumanicus | wc -l)
      EUROPAEUS=$(grep ^Er missingness_$NUM_SAMPLES.txt | cut -f 1 | grep -f - popmap.txt | grep europaeus  | wc -l)
      OTHERS=$(grep ^Er missingness_$NUM_SAMPLES.txt | cut -f 1 | grep -f - popmap.txt | grep -vP "roumanicus|europaeus" | wc -l)
      SITES[0]=$(grep "Sites missing 0 genotypes:" missingness_$NUM_SAMPLES.txt | gawk '{print $5}')
      SITES[1]=$(grep "Sites missing 1 genotypes:" missingness_$NUM_SAMPLES.txt | gawk '{print $5}')
      SITES[2]=$(grep "Sites missing 2 genotypes:" missingness_$NUM_SAMPLES.txt | gawk '{print $5}')
      echo -e "$EXCLUDE\t$NUM_SAMPLES\t$ROUMANICUS\t$EUROPAEUS\t$OTHERS\t$CURRENT_MISSING\t${SITES[2]}\t${SITES[1]}\t${SITES[0]}" >> summary_missingness.txt

      if [ ! -e scores_$NUM_SAMPLES.txt ]; then
         # Vcftools seems to have a bug: the filter --max-missing-count is sometimes more conservative than i should.
         # Thus, I don't trust the number of sites in the 012 files will be the expected from the *.recode.vcf files.
         if [ -e EuropaeusRoumanicus_$NUM_SAMPLES.012.pos ]; then
            NUM_SITES=$(cat EuropaeusRoumanicus_$NUM_SAMPLES.012.pos | wc -l)
         else
            NUM_SITES=0
         fi
         i=0
         while [ $NUM_SITES -le $NUM_SAMPLES ]; do
            vcftools --vcf EuropaeusRoumanicus_$NUM_SAMPLES.recode.vcf --max-missing-count $i --012 --out EuropaeusRoumanicus_$NUM_SAMPLES
            NUM_SITES=$(cat EuropaeusRoumanicus_$NUM_SAMPLES.012.pos | wc -l)
            i=$((i+1))
         done
         sed -i 's/-1/NA/g' EuropaeusRoumanicus_$NUM_SAMPLES.012
         R --save < run_Adegenet.R --args EuropaeusRoumanicus_$NUM_SAMPLES.012 $NUM_SITES EuropaeusRoumanicus_$NUM_SAMPLES.012.indv EuropaeusRoumanicus_$NUM_SAMPLES.012.pos scores_$NUM_SAMPLES.txt
         LC_ALL=C sort -gk 2,2 scores_$NUM_SAMPLES.txt > z1
         mv z1 scores_$NUM_SAMPLES.txt
      fi

      EXCLUDE=$(head -n 1 missingness_$NUM_SAMPLES.txt | cut -f 1)
   done
fi

if [ -e excluded.txt ]; then rm excluded.txt; fi

if [ ! -e coordinates.txt ]; then
   cp $GEODATA ./coordinates.txt
fi

if [ ! -e scores.png ]; then
   if [ ! -e table_scores.txt ]; then
      echo -e "#Num_samples\tSample\tPCA1\PCA2\tLongitude\tLatitude" > table_scores.txt
      for SCORE_FILE in $(ls -1 scores_*.txt); do
         NUM_SAMPLES=${SCORE_FILE:7:2}
         gawk -v NUM_SAMPLES=$NUM_SAMPLES '(FILENAME == "coordinates.txt"){
            LONGITUDE[$1] = $2
            LATITUDE[$1] = $3
         }(FILENAME ~ /^scores_/){
            printf("%u\t%s\t%.6f\t%.6f\t%.2f\t%.2f\n", NUM_SAMPLES, $1, $2, $3, LONGITUDE[$1], LATITUDE[$1])
         }' coordinates.txt $SCORE_FILE >> table_scores.txt
      done
   fi
   R --save < plot_scores.R
fi

# The plot 'scores.png' shows how the scores on the first principal component, which separate E. roumanicus
# from E. europaeus, vary depending on the number of samples included in the computation. In the plot, I have
# normalized the scores, as if 0 means the most extreme E. roumanicus, and 1 is the most extreme (or pure)
# E. europaeus. Note that the colour indicates longitude (east-west geographic coordinate). On the right, with
# a high number of samples, a low number of sites makes the scores unstable. With 36 or less individuals,
# the genetic positions are quite stable, and reliable. Using 36 as the maximum number of samples that provide
# reliable genetic positions, we are using 897 sites without any missing genotype. And the subset of admixed
# individuals would be: Er36_SK24, Er50_R3, Er55_AU7, Er37_SK27, Er74_SP16, and Er72_SIE1B.
