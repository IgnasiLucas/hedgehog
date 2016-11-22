#!/bin/bash

#				2016-11-18
#				----------
#
# I still would like to inform freebayes of the distribution of individuals
# among samples, so that the genotype likelihoods are more accurate. I don't
# think it will make a big difference, but I think it is worth testing.

DATADIR=../../data

if [ ! -e populations.txt ]; then
   echo "Ask Kristýna for the populations"
   exit
fi

# Preliminary explortaion of the qualities and types of variation in the
# megapool.vcf (42GB), shows that a large amount of variants are complex
# combinations of SNPs, insertions, deletions, MNPs, etc. Vcftools does
# not seem to be able to filter those out (I can only see a filter for indels
# based on the lengths of the alleles). Thus, I use grep below to keep only
# snps and indels. Then, I filter by quality.
if [ ! -e snp_indels.vcf.gz ]; then
   if [ ! -e megapool.vcf ]; then
      freebayes -f $DATADIR/reference.fa \
         --vcf megapool.vcf \
         --populations populations.txt \
         --min-mapping-quality 5 \
         --min-base-quality 5 \
         --read-max-mismatch-fraction 0.15 \
         --read-indel-limit 5 \
         --genotype-qualities \
         $DATADIR/*_sortednuc.bam
   fi
#  Before I erase the megapool.vcf file, I will extract some descriptive
#  information that justifies its reduction: the histogram of qualities and
#  the list of different types of variants:
   if [ ! -e megapool_summary.txt ]; then
      gawk 'BEGIN{
         MAX_QUAL=100
      }(/^[^#]/){
         QUAL = sprintf(".0f", $6)
         if (QUAL > MAX_QUAL) MAX_QUAL = QUAL
         QUAL_FREQ[QUAL]++
         split($8,INFO,/;|=/)
         TYPE_FIELD = 0
         for (FIELD in INFO) {
            if (INFO[FIELD] == "TYPE") TYPE_FIELD = FIELD + 1
         }
         TYPE_FREQ[INFO[TYPE_FIELD]]++
      }END{
         print "Histogram of qualities:"
         for (q = 0; q <= MAX_QUAL; q++) {
            if (q in QUAL_FREQ) {
               print q "\t" QUAL_FREQ[q]
            }
         }
         print
         print "Frequencies of types of variants:"
         for (type in TYPE_FREQ) {
            print type "\t" TYPE_FREQ[type]
         }
      }' megapool.vcf > megapool_summary.txt
   fi
   grep -P "^#|TYPE=snp;|TYPE=ins;|TYPE=del;" megapool.vcf > snp_indels.vcf
   vcftools --vcf snp_indels.vcf --out snp_indels_q40 --recode --recode-INFO-all --minQ 40
   mv snp_indels_q40.recode.vcf snp_indels_q40.vcf
   bgzip -c snp_indels_q40.vcf > snp_indels_q40.vcf.gz
   tabix -p vcf snp_indels_q40.vcf.gz
#   rm megapool.vcf
#   rm snp_indels.vcf
#   rm snp_indels_q40.vcf
fi

# Now, I want to see at HWE in the two populations where there are enough
# individuals.

if [ ! -e pop1.hwe ]; then
   if [ ! -e population1.txt ]; then
      gawk '($2 == 1){print $1}' populations.txt > population1.txt
   fi
   vcftools --vcf megapool.vcf --hardy --max-missing 1.0 --out pop1 --keep population1.txt
fi

if [ ! -e pop2.hwe ]; then
   if [ ! -e population2.txt ]; then
      gawk '($2 == 2){print $1}' populations.txt > population2.txt
   fi
   vcftools --vcf megapool.vcf --hardy --max-missing 1.0 --out pop2 --keep population2.txt
fi
