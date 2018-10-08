#!/bin/bash
#
#				2018-09-26
#				==========
#
# I want to re-do the filtering of the vcf file. We need to insert the binary
# presence flag, and then make sure that we select for analysis sites where at least
# four individuals from E. concolor have data.
#
# After giving it some thought, it becomes clear that the Admixture analysis and
# the abba/baba or D tests have different requirements: Admixture estimates allele
# frequencies, and benefits from large numbers of individuals per population. The
# abba/baba test does not seem so sensitive to allele frequencies: even though current
# implementation uses allele frequencies in the sample, it was originally designed to
# work even with single individuals, and it gains power from the number of SNPs used,
# more than from the number of individuals, since it is estimating a proportion of
# SNPs in the genome.
#
# I will produce two filterings, one optimized for Admixture, and one for the abba/baba.

GZVCF=/data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz
ADD_FLAG_DIR=../../bin
MAX_MISSING=63

if [ ! -e NotErinaceus.txt ]; then
   # The Atelerix and Hemiechinus individuals cannot take part in the Admixture analysis,
   # and should be excluded.
   echo Er65_IS25 >  NotErinaceus.txt
   echo Er73_SNG1 >> NotErinaceus.txt
fi

# The binary flags that indicate presence of alleles in individuals are only meaningful
# as long as the individuals keep their order in the vcf file. Thus, any exclusion
# of individuals should be made before the binary flags are added. Thus, I will have to
# add the flags twice: first to the original vcf file, to be able to filter sites by
# composition, and then again after excluding the individuals with too much missing data.
# Note that I need to determine the amount of missing data on a vcf file already filtered
# by composition. Otherwise, there is over 1.4 million sites and most individuals have
# very high levels of missing data.

if [ ! -e Erinaceus_flagged.vcf ]; then
   vcftools --gzvcf $GZVCF \
            --stdout \
            --remove NotErinaceus.txt \
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
   gawk -f $ADD_FLAG_DIR/add_flag.awk > Erinaceus_flagged.vcf
fi

if [ ! -e Erinaceus_r10e10c4.vcf ]; then
   gawk -f filtervcf.awk -v ROU=10 -v EUR=10 -v CON=4 popmap.txt Erinaceus_flagged.vcf > Erinaceus_r10e10c4.vcf
fi

if [ ! -e missingness.txt ]; then
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
   }' Erinaceus_r10e10c4.vcf | sort -nrk 2,2 > missingness.txt
fi

if [ ! -e missing_$MAX_MISSING.txt ]; then
   gawk -v MAX=$MAX_MISSING.txt '($2 > MAX/100){print $1}' missingness.txt > missing_$MAX_MISSING.txt
fi

if [ ! -e Erinaceus_MaxMissing$MAX_MISSING\_flagged.vcf ]; then
   vcftools --gzvcf $GZVCF \
            --stdout \
            --remove NotErinaceus.txt \
            --remove missing_$MAX_MISSING.txt \
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
   gawk -f $ADD_FLAG_DIR/add_flag.awk > Erinaceus_MaxMissing$MAX_MISSING\_flagged.vcf
fi

if [ ! -e ErinMaxMiss$MAX_MISSING\_r10e10c4.vcf ]; then
   gawk -f filtervcf.awk -v ROU=10 -v EUR=10 -v CON=4 popmap.txt Erinaceus_MaxMissing$MAX_MISSING\_flagged.vcf > ErinMaxMiss$MAX_MISSING\_r10e10c4.vcf
fi

if [ ! -e summary.txt ]; then
   for vcf in `ls -1 *.vcf`; do
      echo -e "file\tNumber of sites" > summary.txt
      echo -e $vcf "\t" `grep -P -v '^#' $vcf | wc -l` >> z1
   done
   sort -nrk 2,2 z1 >> summary.txt
   rm z1
   echo -e "\n# Distributions of number of genotyped individuals per site from each population." >> summary.txt
   for vcf in `ls -1 *.vcf`; do
      gawk '(FILENAME ~ /popmap/){
         POP[$1] = $2
         MAX[$2]++
      }(FILENAME ~ /\.vcf/){
         if ($0 ~ /^#CHROM/) {
            for (i = 10; i <= NF; i++) {
               POPi[i] = POP[$i]
            }
         }
         if ($0 ~ /^[^#]/) {
            delete(N)
            for (i = 10; i <= NF; i++) {
               if ($i !~ /\.\/\./) {
                  N[POPi[i]]++
               }
            }
            F["roumanicus"][N["roumanicus"] + 0]++
            F["europaeus"][N["europaeus"] + 0]++
            F["concolor"][N["concolor"] + 0]++
            F["Hemiechinus"][N["Hemiechinus"] + 0]++
            F["wrong"][N["wrong"] + 0]++
            F["hybrid"][N["hybrid"] + 0]++
            F["Atelerix"][N["Atelerix"] + 0]++
         }
      }END{
         print "\n# " FILENAME
         print "# Number\troumanicus\teuropaeus\tconcolor\tHemiechinus\tAtelerix\thybrid\twrong"
         for (n = 0; n <= 24; n++) {
            print n "\t" F["roumanicus"][n]+0 "\t" F["europaeus"][n]+0 "\t" F["concolor"][n]+0 "\t" F["Hemiechinus"][n]+0 "\t" F["Atelerix"][n]+0 "\t" F["hybrid"][n]+0 "\t" F["worng"][n]+0
         }
      }' popmap.txt $vcf >> summary.txt
   done
fi
