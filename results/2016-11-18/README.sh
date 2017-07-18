#!/bin/bash
#				2016-11-18
#				----------
#
# I still would like to inform freebayes of the distribution of individuals
# among samples, so that the genotype likelihoods are more accurate. I don't
# think it will make a big difference, but I think it is worth testing.

DATADIR=../../data

if [ ! -e $DATADIR/populations.txt ]; then
   echo "Ask Kristýna for the populations"
   exit
fi

# Preliminary explortaion of the qualities and types of variation in the
# megapool.vcf (42GB), shows that a large amount of variants include more
# than one alternative alleles, each of potentially many different kinds.
# I considered retaining only biallelic variants. However, given that there
# are relatively few polymorphic sites, and because the divergence between
# populations (actually, species) may allow for multiple alleles in several
# loci, I am tempted to include at least some tri-allelic sites. Here, I
# will use biallelic snps and indels. For comparison, I may run the whole
# analysis in a separate folder with triallelic variants.

if [ ! -e megapool_summary.txt ]; then
   if [ ! -e megapool.vcf ]; then
      freebayes -f $DATADIR/reference.fa \
         --vcf megapool.vcf \
         --populations $DATADIR/populations.txt \
         --theta 0.01 \
         --standard-filters \
         --max-coverage 900 \
         --genotype-qualities \
         $DATADIR/*.bam

#         --min-mapping-quality 5 \
#         --min-base-quality 5 \
#         --read-max-mismatch-fraction 0.15 \
#         --read-indel-limit 5 \
   fi

#  Before I erase the megapool.vcf file, I will extract some descriptive
#  information that justifies its reduction: the histogram of qualities and
#  the list of different types of variants:
   gawk 'BEGIN{
      MAX_QUAL=100
   }(/^[^#]/){
      QUAL = sprintf("%.0f", $6)
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
      print ""
      print "Frequencies of types of variants:"
      for (type in TYPE_FREQ) {
         print type "\t" TYPE_FREQ[type]
      }
   }' megapool.vcf > megapool_summary.txt
fi

# These are the types of variants making up the 99% of the variable sites:
#
# ------------------------------------------------
#  Type                  Number   Accumulated (%)
# ------------------------------------------------
#  snp                  2332548           92.4759
#  complex                58484           94.7945
#  mnp                    48885           96.7326
#  snp,complex            24729           97.7130
#  snp,snp                15998           98.3472
#  snp,mnp                13579           98.8856
#  complex,snp            13535           99.4222
#  mnp,snp                 7280           99.7108
#  complex,complex         1329           99.7635
# ------------------------------------------------
#
# It is unclear why there are no indels. It could be that the standard filters
# applied in the last run remove indels, although that would be unexpected.

if [ ! -e snp_indels_q40.vcf.gz ]; then
   if [ -e megapool.vcf ]; then
      grep -P "^#|TYPE=snp[;\t]|TYPE=ins[;\t]|TYPE=del[;\t]" megapool.vcf > snp_indels.vcf
      vcftools --vcf snp_indels.vcf --minQ 40 --out snp_indels_q40 --recode --recode-INFO-all
      mv snp_indels_q40.recode.vcf snp_indels_q40.vcf
      bgzip -c snp_indels_q40.vcf > snp_indels_q40.vcf.gz
      tabix -p vcf snp_indels_q40.vcf.gz
#     rm megapool.vcf
#     rm snp_indels.vcf
#     rm snp_indels_q40.vcf
   else
      echo "megapool.vcf is missing. I suggest you remove megapool_summary.txt and run README.sh again."
      exit
   fi
fi

# Now, I want to see at HWE in the two populations where there are enough
# individuals.

if [ ! -e summary_snp_pop1.hwe ]; then
   if [ ! -e population1.txt ]; then
      gawk '($2 == 1){print $1}' $DATADIR/populations.txt > population1.txt
   fi
   vcftools --gzvcf snp_indels_q40.vcf.gz --max-missing 1.0 --keep population1.txt --remove-indels --out snp1 --hardy
   echo -e "OBS(HOM1/HET/HOM2)\tP_HWE\tP_HET_DEFICIT\tP_HET_EXCESS\tFREQUENCY" > summary_snp_pop1.hwe
   gawk '(NR > 1){
      FREQ[$3 "\t" $6 "\t " $7 "\t" $8]++
   }END{
      for (i in FREQ) print i "\t" FREQ[i]
   }' snp1.hwe | sort -nrk 5 >> summary_snp_pop1.hwe
   rm snp1.hwe
fi

if [ ! -e summary_snp_pop2.hwe ]; then
   if [ ! -e population2.txt ]; then
      gawk '($2 == 2){print $1}' $DATADIR/populations.txt > population2.txt
   fi
   vcftools --gzvcf snp_indels_q40.vcf.gz --max-missing 1.0 --keep population2.txt --remove-indels --out snp2 --hardy
   echo -e "OBS(HOM1/HET/HOM2)\tP_HWE\tP_HET_DEFICIT\tP_HET_EXCESS\tFREQUENCY" > summary_snp_pop2.hwe
   gawk '(NR > 1){
      FREQ[$3 "\t" $6 "\t" $7 "\t" $8]++
   }END{
      for (i in FREQ) print i "\t" FREQ[i]
   }' snp2.hwe | sort -nrk 5 >> summary_snp_pop2.hwe
   rm snp2.hwe
fi

if [ ! -e snp_indels_q40.het ]; then
   vcftools --gzvcf snp_indels_q40.vcf.gz --het --out snp_indels_q40
   gawk -v OFS="\t" '(FILENAME ~ /populations.txt/){
      POP[$1] = $2
   }(FILENAME ~ /snp_indels_q40\.het/){
      print $0 "\t" POP[$1]
   }' $DATADIR/populations.txt snp_indels_q40.het > z1
   mv z1 snp_indels_q40.het
fi

# The following should have sent standard errors to individual files.
for i in 1 2; do
   if [ ! -e pop$i.het ]; then
      vcftools --gzvcf snp_indels_q40.vcf.gz --out pop$i --keep population$i.txt --het &
   fi
   if [ ! -e pop$i.relatedness ]; then
      vcftools --gzvcf snp_indels_q40.vcf.gz --out pop$i --keep population$i.txt --relatedness &
   fi
   if [ ! -e pop$i.idepth ]; then
      vcftools --gzvcf snp_indels_q40.vcf.gz --out pop$i --keep population$i.txt --depth &
   fi
#   if [ ! -e pop$i.geno.ld ]; then
#      vcftools --gzvcf snp_indels_q40.vcf.gz --out pop$i --keep population$i.txt --geno-r2 &
#   fi
#   if [ ! -e pop$i.interchrom.geno.ld ]; then
#      vcftools --gzvcf snp_indels_q40.vcf.gz --out pop$i --keep population$i.txt --interchrom-geno-r2 &
#   fi
done
#wait

if [ ! -e pop1_pop2.weir.fst ]; then
   gawk '($2 ~ /3|4|5|6/){print $1}' $DATADIR/populations.txt > populations3to6.txt
   vcftools --gzvcf snp_indels_q40.vcf.gz --out pop1_pop2 --remove populations3to6.txt --weir-fst-pop population1.txt &
fi
wait

# CONCLUSIONS
# -----------
#
# From the HWE statistics, it is clear that most of the 'variable' sites are not so
# variable within a specific population. Among sites with complete data in population
# 1 (romanicus, 24 individuals), ~50% of the 1276 SNPs (identified among all populations)
# are monomorfic for either allele. In population 2 (20 individuals), among the 376 SNPs
# with complete data, ~47% are monomorfic for the either allele. In addition, we observe
# a deficit of heterozygotes in more than 17% of sites in both populations (at 0.05
# significance level).
#
