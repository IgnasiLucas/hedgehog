#!/bin/bash
#
#				2020-04-07
#				==========
#
# One of the open issues in this research is the origin of the negative correlation
# between genetic diversity and divergence in the genome of E. roumanicus, when
# divergence is measured against E. europaeus (2018-10-29). Under neutrality, a
# positive correlation is expected (Kimura, 1983). That is the case between
# D. melanogaster and D. simulans (Begun et al. 2017, https://dx.doi.org/10.1371%2Fjournal.pbio.0050310).
# Also among 40 strains of S. cerevisiae (Cutter & Moses, 2011, https://doi.org/10.1093/molbev/msq356).
# Interestingly, Cutter and Moses (2011) distinguish neutral- from replacement-site
# divergence and polymorphism. They do observe the neutral expectation of positive
# correlation between divergence and polymorphism, but only in non-neutral,
# replacement sites. Among sites assumed neutral the correlation is very weak
# at best, and of different sign at different genomic scales (see fig. 2 in
# Cutter and Moses, 2011).
#
# This is to stress that the lack of correlation between polymorphism and
# divergence in E. europaeus also departs from the neutral expectation and
# needs to be explained.
#
# Our co-author Pavel Hulva suggests including E. concolor in this analysis. The
# reason is that we do not detect significant introgression between E. roumanicus
# and E. concolor. Thus, if a negative correlation persists in E. roumanicus between
# polimorphism and divergence, when divergence is measured against E. concolor, then
# we could no longer attribute it to introgression.
#
# We still need to take into account introgression from E. europaeus into E. roumanicus.
# Because E. concolor is closer to E. roumanicus than to E. europaeus, introgressed
# variation from E. europaeus into E. roumanicus can also inflate divergence
# between E. roumanicus and E. concolor. I will just repeat the whole analysis
# integrating E. concolor.
#
#------------------------------------------------------
# SETTINGS
#------------------------------------------------------

FREQ2_AWK=../../bin/freq2.awk

# Different vcf files are optimized for different statistics:

VCF1=../2018-09-26/ErinaceusAndHemiechinus_r1e1c1h1.vcf
VCF2=../2018-09-26/ErinMaxMiss63_r10e10c4.vcf

# This is the folder where Simon H. Martin's genomic_general repository is cloned.
# S.H. Martin added some support for python 3. But some scripts still need python 2
# and numpy.

GENOMICS_GENERAL=~/bin/genomics_general
POPDATA=../../data/populations.txt

#-------------------------------------------------------
# POPULATION DATA
#-------------------------------------------------------

# I will exclude both the hybrid and the other admixed individual, in order to work only
# with the historical signal of introgression. It will help interpretation. The hybrid
# individual is already identified as such in the popmap file. Below, I assign Er55_AU7
# to a population called 'admixed'. That's enough for them to be ignored when running
# either popgenWindows.py or ABBABABAwindows.py.

if [ ! -e creh1.geno.gz ]; then
   python $GENOMICS_GENERAL/VCF_processing/parseVCF.py -i $VCF1 | gzip > creh1.geno.gz
fi


if [ ! -e erin63.geno.gz ]; then
   python $GENOMICS_GENERAL/VCF_processing/parseVCF.py -i $VCF2 | gzip > erin63.geno.gz
fi

if [ ! -e creh1.popmap.txt ] || [ ! -e erin63.popmap.txt ]; then
   if [ ! -e popmap.txt ]; then
      gawk 'BEGIN{
         POPNAME[1] = "roumanicus"
         POPNAME[2] = "europaeus"
         POPNAME[3] = "concolor"
         POPNAME[4] = "hybrid"
         POPNAME[5] = "Hemiechinus"
         POPNAME[6] = "Atelerix"
      }(/^[^#]/){
         if ($1 != "Er55_AU7") {
            print $1 "\t" POPNAME[$2]
         } else {
            print $1 "\tadmixed"
         }
      }' $POPDATA | sort -k 2,2 > popmap.txt
   fi

   if [ ! -e creh1.popmap.txt ]; then
      gunzip -c creh1.geno.gz | head -n 1 | sed 's/\t/\n/g' | tail -n +3 > good_samples.txt
      grep -F -f good_samples.txt popmap.txt > creh1.popmap.txt
      rm good_samples.txt
   fi
   if [ ! -e erin63.popmap.txt ]; then
      gunzip -c erin63.geno.gz | head -n 1 | sed 's/\t/\n/g' | tail -n +3 > good_samples.txt
      grep -F -f good_samples.txt popmap.txt > erin63.popmap.txt
      rm good_samples.txt
      rm popmap.txt
   fi
fi

#---------------------------------------------------------
#                WHOLE GENOME STATISTICS
#---------------------------------------------------------
#
# I want to estimate divergence, polymorphism and introgression in the same non-overlapping
# windows. Introgression statistics are more demanding on the number of SNPs than
# population genetic ones. By obtaining first introgression statistics, the windows are set
# to accommodate their requirements, and then re-used in the estimation of population genetics
# statistics.

if [ ! -e nonoverlap.abbababa.csv ]; then
   python $GENOMICS_GENERAL/ABBABABAwindows.py -g creh1.geno.gz -f phased -o nonoverlap.abbababa.csv --windType sites \
   -w 110 --overlap 0 -m 100 --popsFile creh1.popmap.txt -P1 concolor -P2 roumanicus -P3 europaeus -O Hemiechinus --writeFailedWindows -T 20
fi
if [ ! -e nonoverlap.PopGenStats.csv ]; then
   if [ ! -e nonoverlap.windows.txt ]; then
      cut -f 1,2,3 -d ',' --output-delimiter=$'\t' nonoverlap.abbababa.csv | tail -n +2 > nonoverlap.windows.txt
   fi
   python $GENOMICS_GENERAL/popgenWindows.py --windType predefined --windCoords nonoverlap.windows.txt -f phased \
   -g erin63.geno.gz -o nonoverlap.PopGenStats.csv --popsFile erin63.popmap.txt -p roumanicus -p europaeus -p concolor -T 20 --writeFailedWindows
fi

#-----------------------------------------------------------
# GENES AND INTERGENIC REGIONS
#-----------------------------------------------------------
#
# Genes and intergenic regions are disjoint windows, but I can't separate them using regular windows, and I want
# to estimate parameters separately for the two kinds of regions. Thus, below I split the geno files. Warning:
# I assume that whatever does not overlap a gene is intergenic.
#
# I already estimated the genome-wide introgression statistics for genic and intergenic regions
# in 2018-10-29. I do not need them here. Thus, I only split the erin63.geno.gz file.

if [ ! -e ../../data/annotation.gff3.gz ]; then
   wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Erinaceus_europaeus/GFF/ref_EriEur2.0_top_level.gff3.gz
   mv ref_EriEur2.0_top_level.gff3.gz ../../data/annotation.gff3.gz
fi

if [ ! -e genes.bed ]; then
   gunzip -c ../../data/annotation.gff3.gz | grep -P "^#|\tgene\t" | bedtools merge -i - > genes.bed
fi

if [ ! -e genes.erin63.geno.gz ]; then
   gunzip -c erin63.geno.gz | head -n 1 > genes.erin63.geno
   if [ ! -e genic_filter.erin63.txt ]; then
      if [ ! -e snps.erin63.bed ]; then
         gunzip -c erin63.geno.gz | gawk '(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' > snps.erin63.bed
      fi
      bedtools intersect -a snps.erin63.bed -b genes.bed -wa | gawk '{print $1 "\t" $3 "\t"}' > genic_filter.erin63.txt
   fi
   gunzip -c erin63.geno.gz | grep -F -f genic_filter.erin63.txt >> genes.erin63.geno
   gzip genes.erin63.geno
fi

if [ ! -e inter.erin63.geno.gz ]; then
   # No need to add a header. It gets added latter from erin63.geno.gz, because the filter is reversed now.
   #gunzip -c erin63.geno.gz | head -n 1 > inter.erin63.geno
   if [ ! -e genic_filter.erin63.txt ]; then
      if [ ! -e snps.erin63.bed ]; then
         gunzip -c erin63.geno.gz | gawk '(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' > snps.erin63.bed
      fi
      bedtools intersect -a snps.erin63.bed -b genes.bed -wa | gawk '{print $1 "\t" $3 "\t"}' > genic_filter.erin63.txt
   fi
   gunzip -c erin63.geno.gz | grep -F -v -f genic_filter.erin63.txt >> inter.erin63.geno
   gzip inter.erin63.geno
   rm genic_filter.erin63.txt snps.erin63.bed
fi

if [ ! -e genes.PopGenStats.csv ]; then
   python $GENOMICS_GENERAL/popgenWindows.py --windType sites -w 110 --overlap 0 -m 100 -f phased \
      -g genes.erin63.geno.gz -o genes.PopGenStats.csv --popsFile erin63.popmap.txt -p roumanicus -p europaeus -p concolor -T 20 --writeFailedWindows
fi

if [ ! -e inter.PopGenStats.csv ]; then
   python $GENOMICS_GENERAL/popgenWindows.py --windType sites -w 110 --overlap 0 -m 100 -f phased \
      -g inter.erin63.geno.gz -o inter.PopGenStats.csv --popsFile erin63.popmap.txt -p roumanicus -p europaeus -p concolor -T 20 --writeFailedWindows
fi

#------------------------------------------------------------------
#                             PLOTS
#------------------------------------------------------------------

if [ ! -e HKA.pdf ]; then
   R -q --no-save -e "rmarkdown::render('HKA.Rmd', output_file='HKA.pdf')"
fi
