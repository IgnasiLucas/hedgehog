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
# -------------------------------------
# Type               Number  Proportion
# -------------------------------------
# snp              19752403   0.754746
# complex           1625297   0.816849
# ins               1068231   0.857666
# del               1032943   0.897135
# mnp                776732   0.926814
# complex,snp        289130   0.937862
# snp,complex        274333   0.948345
# snp,snp            225525   0.956962
# mnp,snp            127073   0.961817
# snp,mnp            115994   0.96625
# del,snp             97237   0.969965
# complex,complex     70028   0.972641
# ins,ins             62622   0.975034
# ins,snp             61884   0.977398
# snp,del             52852   0.979418
# del,del             37097   0.980835
# snp,ins             32865   0.982091
# del,ins             30788   0.983267
# del,complex         28069   0.98434
# complex,del         22828   0.985212
# complex,mnp         22221   0.986061
# ins,del             20750   0.986854
# ins,complex         17124   0.987508
# snp,complex,complex 17060   0.98816
# complex,ins         15651   0.988758
# snp,complex,snp     14898   0.989328
# complex,snp,snp     13920   0.989859
# mnp,mnp             11245   0.990289
# -------------------------------------
#
# SNPs, insertions and deletions consitute less than 84% of the variable sites.

if [ ! -e snp_indels_q40.vcf.gz ]; then
   if [ -e megapool.vcf ]; then
      grep -P "^#|TYPE=snp;|TYPE=ins;|TYPE=del;" megapool.vcf > snp_indels.vcf
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

if [ ! -e summary_indel_pop1.hwe ]; then
   if [ ! -e population1.txt ]; then
      gawk '($2 == 1){print $1}' $DATADIR/populations.txt > population1.txt
   fi
   vcftools --gzvcf snp_indels_q40.vcf.gz --max-missing 1.0 --keep population1.txt --keep-only-indels --out indel1 --hardy
   echo -e "OBS(HOM1/HET/HOM2)\tP_HWE\tP_HET_DEFICIT\tP_HET_EXCESS\tFREQUENCY" > summary_indel_pop1.hwe
   gawk '(NR > 1){
      FREQ[$3 "\t" $6 "\t" $7 "\t" $8]++
   }END{
      for (i in FREQ) print i "\t" FREQ[i]
   }' indel1.hwe | sort -nrk 5 >> summary_indel_pop1.hwe
   rm indel1.hwe
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

if [ ! -e summary_indel_pop2.hwe ]; then
   if [ ! -e population2.txt ]; then
      gawk '($2 == 2){print $1}' $DATADIR/populations.txt > population2.txt
   fi
   vcftools --gzvcf snp_indels_q40.vcf.gz --max-missing 1.0 --keep population2.txt --keep-only-indels --out indel2 --hardy
   echo -e "OBS(HOM1/HET/HOM2)\tP_HWE\tP_HET_DEFICIT\tP_HET_EXCESS\tFREQUENCY" > summary_indel_pop2.hwe
   gawk '(NR > 1){
      FREQ[$3 "\t" $6 "\t" $7 "\t" $8]++
   }END{
      for (i in FREQ) print i "\t" FREQ[i]
   }' indel2.hwe | sort -nrk 5 >> summary_indel_pop2.hwe
   rm indel2.hwe
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
# 1 (22 individuals), 70% of the 23585 SNPs identified among all populations are
# monomorfic for the reference allele. Among indels 55% are also fixed for the reference
# allele in population 1. In population 2 (20 individuals), among the 1345 SNPs with
# complete data, 54% are monomorfic for the reference allele; and 24% of the 206 indels
# are also fixed for the reference. In addition, we observe a deficit of heterozygotes
# in more than 15% of sites in both populations (at 0.05 significance level).
#
# Curiously, the estimates of heterozygosity and inbreeding coefficient obtained for
# all individuals together suggest that populations 1, and 2 have an excess of heterozygote
# sites (negative inbreeding coefficient of around -0.1), and the same is true for the
# individual of population 3 (F=-0.3) while individuals from populations 5 and 6 seem
# to be highly inbred (F=0.48, F=0.36, respectively). However, these values are based
# on only 440 sites. Files pop1.het and pop2.het estimate heterozygosity and F for populations
# 1 and 2 separately. Then, 7793 and 695 sites are used respectively, and the median
# inbreeding coefficients become positive: 0.34 for population 1 and 0.02 for population 2.
#
