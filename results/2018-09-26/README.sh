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
INDLIST=/data/kristyna/hedgehog/results_2018/05-06-2018/r10e10c4/popmap_erinaceus

if [ ! -e IndList.txt ]; then
   # This file is just the list of samples to be used, after removing some with low
   # quality data and excluding the two samples from diferen genera (Atelerix and Hemiechinus).
   cut -f 1 $INDLIST > IndList.txt
fi
for rou in 10 15; do
   for eur in 10 15; do
      for con in 4 5; do
         if [ ! -e r${rou}e${eur}c${con}.vcf ]; then
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
            gawk -f filtervcf.awk -v ROU=$rou -v EUR=$eur -v CON=$con popmap.txt flagged.vcf > r${rou}e${eur}c${con}.vcf &
         fi
      done
   done
done
wait

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
#            N["roumanicus"] = 0; N["europaeus"] = 0; N["concolor"] = 0; N["Hemiechinus"] = 0; N["wrong"] = 0
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
