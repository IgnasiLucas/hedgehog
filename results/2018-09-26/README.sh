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

if [ ! -e AtelerixAndWrong.txt ]; then
   # The Atelerix and Hemiechinus individuals cannot take part in the Admixture analysis,
   # and should be excluded. Also, Er59_FR1 and Er68_PRT1B seem to have problems. However,
   # the Hemiechinus individual is used as an outgroup in abba/baba test, and needs to be
   # in the vcf file used for that analysis.
#   echo Er65_IS25  >  NotErinaceus.txt    # Hemiechinus
   echo Er73_SNG1  >> AtelerixAndWrong.txt    # Atelerix
   echo Er59_FR1   >> AtelerixAndWrong.txt    # wrong
   echo Er68_PRT1B >> AtelerixAndWrong.txt    # wrong
fi

# The binary flags that indicate presence of alleles in individuals are only meaningful
# as long as the individuals keep their order in the vcf file. Thus, any exclusion
# of individuals should be made before the binary flags are added. Thus, I will have to
# add the flags twice: first to the original vcf file, to be able to filter sites by
# composition, and then again after excluding the individuals with too much missing data.
# Note that I need to determine the amount of missing data on a vcf file already filtered
# by composition. Otherwise, there is over 1.4 million sites and most individuals have
# very high levels of missing data.

if [ ! -e ErinaceusAndHemiechinus_flagged.vcf ]; then
   vcftools --gzvcf $GZVCF \
            --stdout \
            --remove AtelerixAndWrong.txt \
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
   gawk -f $ADD_FLAG_DIR/add_flag.awk > ErinaceusAndHemiechinus_flagged.vcf
fi

if [ ! -e ErinaceusAndHemiechinus_r10e10c4.vcf ]; then
   # This is vcf file is just an intermediate file to determine the proportion of missing
   # sites per individual. Not to be used in subsequent analyses. It can be erased.
   gawk -f filtervcf.awk -v ROU=10 -v EUR=10 -v CON=4 popmap.txt ErinaceusAndHemiechinus_flagged.vcf > ErinaceusAndHemiechinus_r10e10c4.vcf
fi

if [ ! -e ErinaceusAndHemiechinus_r1e1c1h1.vcf ]; then
   # This is the vcf file meant to be used for the abba/baba test, including only sites
   # with at least 1 Hemiechinus, 1 roumanicus, 1 europaeus, and 1 concolor.
   gawk -f filtervcf.awk -v HEM=1 -v ROU=1 -v EUR=1 CON=1 popmap.txt ErinaceusAndHemiechinus_flagged.vcf > ErinaceusAndHemiechinus_r1e1c1h1.vcf
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
   }' ErinaceusAndHemiechinus_r10e10c4.vcf | sort -nrk 2,2 > missingness.txt
fi

if [ ! -e missing_$MAX_MISSING.txt ]; then
   gawk -v MAX=$MAX_MISSING.txt '($2 > MAX/100){print $1}' missingness.txt > missing_$MAX_MISSING.txt
fi

if [ ! -e NotErinaceus.txt ]; then
   cp AtelerixAndWrong.txt NotErinaceus.txt
   echo Er65_IS25      >>  NotErinaceus.txt
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
