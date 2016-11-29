#!/bin/bash
#				2016-11-24
#				----------
#
# I am looking at the variation discovered in the alignments of the reads to
# the reference genome with freebayes. The strategy has been to generate first
# a vcf file with all potentially interesting sites, being permissive with the
# quality filters. We can always filter at later steps of the analysis. In the
# folder 2016-11-18, I created the megapool.vcf and selected biallelic sites
# with either a SNP or an indel, and with a minimum quality of 40. Here, I will
# add complex variants and multiple nucleotide variants, and also triallelic
# sites. The rationale is that there are relatively few variants shared among
# populations, because they are different species. Given the distance between
# populations, I expect triallelic variants to contribute relatively more to
# the amount of shared variants than to the total number of variants. That is,
# the previous filtering of triallelic variants could have removed proportionally
# more shared variants than population-specific ones.
#

DATADIR=../../data
LASTDIR=`pwd | sed 's/2016-11-24/2016-11-18/'`

if [ ! -e $DATADIR/populations.txt ]; then
   echo "Ask Kristýna for the populations"
   exit
fi

# If by the time this runs the 2016-11-18/megapool.vcf has not been removed, I
# will use that. Otherwise, I will have to create it again. Now, I will name it
# something else:

if [ ! -e bitriallelic40.vcf.gz ]; then
   if [ ! -e hedgehog.vcf ]; then
      if [ -e $LASTDIR/megapool.vcf ]; then
         # Note that this is intentionally a hard link, because the target is to be erased.
         ln $LASTDIR/megapool.vcf ./hedgehog.vcf
      else
         freebayes -f $DATADIR/reference.fa \
            --vcf hedgehog.vcf \
            --populations $DATADIR/populations.txt \
            --min-mapping-quality 5 \
            --min-base-quality 5 \
            --read-max-mismatch-fraction 0.15 \
            --read-indel-limit 5 \
            --genotype-qualities \
            $DATADIR/*_sortednuc.bam
      fi
   fi

   # Here, I use gawk to filter the vcf file. This is a little bit risky.
   if [ ! -e bitriallelic40.vcf ]; then
      gawk -f filter.awk hedgehog.vcf > bitriallelic40.vcf
   fi
   bgzip bitriallelic40.vcf
   tabix -p vcf bitriallelic40.vcf.gz
   #rm hedgehog.vcf
fi

# I want to make sure that the distribution of qualities are similar among types of
# variants. Note that old bcftools does not read vcf.gz files, but version 1.3.1 does.
if [ ! -e qualities.png ]; then
   ~/bin/bcftools-1.3.1/bcftools view -H bitriallelic40.vcf.gz | gawk '{
      split($8,A,/TYPE=/)
      split(A[2],B,/;/)
      print $6 "\t" B[1]
   }' > z1
   R --no-save < plot_qualities.R
   rm z1
fi

# I will get HWE summaries also here, to compare with the previous analysis, but without
# distinguishing the type of variation, because now we keep many more types. In any case,
# I assume HWE statistics may only be available for biallelic variants.

for i in 1 2; do
   if [ ! -e summary_pop$i.hwe ]; then
      if [ ! -e population$i.txt ]; then
         gawk -v POP=$i '($2 == POP){print $1}' $DATADIR/populations.txt > population$i.txt
      fi
      if [ ! -e pop$i.hwe ]; then
         vcftools --gzvcf bitriallelic40.vcf.gz --max-missing 1.0 --keep population$i.txt --out pop$i --hardy
      fi
      echo -e "OBS(HOM1/HET/HOM2)\tP_HWE\tP_HET_DEFICIT\tP_HET_EXCESS\tFREQUENCY" > summary_pop$i.hwe
      gawk '(NR > 1){
         FREQ[$3 "\t" $6 "\t " $7 "\t" $8]++
      }END{
         for (i in FREQ) print i "\t" FREQ[i]
      }' pop$i.hwe | sort -nrk 5 >> summary_pop$i.hwe
      rm pop$i.hwe
   fi
done


