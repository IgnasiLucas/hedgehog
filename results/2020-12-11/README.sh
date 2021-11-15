#!/bin/bash
#
#				2020-12-11
#				==========
#
# The Jackknife
# =============
#
# In computing D and f statistics, what is the sampling unit? The single site
# configuration across the species. The statistics are based on counts of specific
# configurations, and their proportions. But, because sites are not independent,
# the sampled configurations are not independent identically distributed (i.i.d.)
# variables. Linkage is the reason for dependence among sites: the same configuration
# is likely to appear among neighbour sites because they share the same genealogy,
# up to the last time recombination happened between them.
#
# When using one individual per species and counting BABA and ABBA patterns,
# variation in the counts is due to the random sampling of individuals and to
# the evolutionary variance among sites. When using not one, but several individuals
# per species, instead of observing individual site patterns, we estimate their
# frequencies (per site). The variance in D and f due to random sampling of individuals
# is reduced, even though there is still variation in every site's esimate of
# configuration frequencies, related to the number of individuals sampled per
# population.
#
# What about the among-sites dependence structure? Do several individuals per
# population reduce the dependence induced by linkage? It seems to me that if
# all populations were in linkage equilibrium, genealogies would be independent
# among neighbour sites. However, disequilibrium happens, especially at short
# distances. And admixture is an expected cause of LD. Thus, a block jackknife
# approach seems to be still necessary.
#
# The inclusion of recently admixed individuals in the sample introduces linkage
# disequilibrium along large portions of some chromosomes. Larger blocks may be
# required to re-sample sites in the jackknife procedure.
#

ADDBPF='../../bin/add_flag.awk'
FREQ3='../../bin/freq3.awk'
ADMIXTURE_DIR='/data/kristyna/hedgehog/results_2020/2020-06-12'
POPMAP_DIR='/data/kristyna/hedgehog/results_2020/2020-04-20'
SOURCE_VCF='/data/kristyna/hedgehog/results_2020/2020-02-28/all_merged_2020.vcf'
DATADIR=$(pwd | sed 's/results/data/')
LINEAGE=( 'E_roumanicus_Asian' 'E_europaeus_Iberian' 'E_concolor' \
          'E_amurensis' 'E_roumanicus_Balkan' 'E_europaeus_Apennine' )

# The $LINEAGE array is manually edited to indicate in every position the lineage
# that corresponds with that position in the Admixture Q file, with K=6. Should be
# checked with EH.popmap_check.txt file, below. Because it will be passed to awk,
# I am using position 0 as 1, intentionally. See below.

# =================================================================================
#  PREPARING THE DATA
# =================================================================================
#
# # POPULATION ASSIGNMENT
#   =====================
#
# The original species assignment of the two sample sets are in ../../data. But I need
# to contrast that with the assignments by Admixture, which is the only source of the
# assignment to the two E. roumanicus clades. According to the manuscript, the number
# of groups must be 6, only Erinaceus samples are used in the analysis. I copy Kristyna's
# files, to make sure I keep a copy. The vcf file must be the one that includes also
# the outgroup species. Below I depart from the vcf file produced with freebayes in
# Kristyna's 2020-02-28 folder. It needs filtering.
#
# I could not keep from checking the filtering settings. A histogram of log-site depth
# has a bump at a site depth between 80 and 90, drops rapidly until depth 200, and then
# has a second smaller bump until depth 1000. We have been using a --maxDP 200 setting,
# under the impression that the first, large bump is the bulk of good quality data.
# There are two mistakes here: --maxDP sets a maximum depth for a genotype, not for the
# site; and actually a total depth around 85 among 90 samples is quite poor. The fact
# is that there is a huge amount of sites that need to be removed for their low
# coverage. It is the second bump, with total depth up to 1000 that must have the
# good data. Luckily, my first mistake makes up for the second, and the resulting vcf
# file is dominated by total depths arround 800.
#
# Other filters  are also common to those used by Kristyna in 2020-02-28, the difference
# being the list of samples to keep, because I need the outgroup, Hemiechinus. Note that
# for the purpose of the present analysis, I include hybrid individuals as well.

if [ ! -d $DATADIR ]; then mkdir $DATADIR; fi
if [ ! -e $DATADIR/admixture.K6.Q ]; then
   cp $ADMIXTURE_DIR/erinaceus_all_2020.6.Q $DATADIR/admixture.K6.Q
fi

# admixture.K6.Q shows the proportion of each of the 6 presumed ancestries in the 75
# individuals included in the Admixture analysis. Below, the admixture.names.txt shows
# the names of individuals in the same order.

if [ ! -e $DATADIR/admixture.names.txt ]; then
   cp $ADMIXTURE_DIR/erinaceus_merged_2020_only_good.012.indv $DATADIR/admixture.names.txt
fi

if [ ! -e $DATADIR/popmap.txt ]; then
   cp $POPMAP_DIR/popmap_all_correct $DATADIR/popmap.txt
fi

# Note that the popmap.txt file includes all 75 samples in the Admixture results, and some
# more. It makes sense to limit the present analysis to the samples that we can assign to
# specific populations on the grounds of the Admixture analysis, even though I would miss
# some individuals with low coverage that were excluded from the Admixture. It is better
# to be consistent with Admixture and skip those. The only individual I need to add with
# respect to the Admixture analysis is the outgroup, Hemiechinus auritus, namely Er65_IS25.
# Unless I remove them, there will be two Erinaceus amurensis individuals, which are a sister
# clade of the E. europaeus species, wihtout any expected hybridization with E. roumanicus...

if [ ! -e EH_list.txt ]; then
   cat $DATADIR/admixture.names.txt > EH_list.txt
   echo "Er65_IS25" >> EH_list.txt
fi
if [ ! -e EH.vcf ]; then
   if [ ! -e EH.recode.vcf ]; then
      vcftools --vcf $SOURCE_VCF \
               --keep EH_list.txt \
               --recode \
               --maf 0.0125 \
               --remove-indels \
               --min-alleles 2 \
               --max-alleles 2 \
               --maxDP 100 \
               --minDP 4 \
               --minQ 50 \
               --thin 261 \
               --max-missing 0.5 \
               --recode-INFO-all \
               --out EH 2> vcftools.log
   fi
   # CAUTION! Unless used with --bignum, the added binary presence flags would be wrong.
   # There are 76 individuals, which produces very big numbers, rounded by default without
   # the --bignum flag.
   gawk --bignum -f $ADDBPF EH.recode.vcf > EH.vcf
   #rm EH.recode.vcf
fi

# The original popmap.txt does not have the information, provided by admixture,
# of the lineages within species. Plus, I need to assign the hybrids to whichever
# group contributes the most to their genome.
if [ ! -e EH.popmap.txt ]; then
   paste $DATADIR/admixture.names.txt $DATADIR/admixture.K6.Q > admixture
   gawk '(FILENAME ~ /admixture/){
      ANCESTRY = 0
      GROUP[$1] = 0
      for (i = 1; i <= 6; i++) {
         if ($(i+1) > ANCESTRY) {
            ANCESTRY = $(i+1)
            GROUP[$1] = i
         }
      }
   }((FILENAME ~ /popmap/) && (/^Er/)){
      SPECIES[$1] = $2
   }((FILENAME ~ /popmap/) && (/^#/)){
      NAME[$2] = $4 "_" $5
   }((FILENAME ~ /EH_list.txt/)){
      print $1 "\t" NAME[SPECIES[$1]] "\t" GROUP[$1]
   }' admixture $DATADIR/popmap.txt EH_list.txt | sort -k 2,3 | \
   tee EH.popmap_check.txt | \
   gawk -v LINEAGE="${LINEAGE[*]}" 'BEGIN{
      split(LINEAGE, LINARRAY, " ")
      # Only Hemiechinus expected to lack group assignment.
      LINARRAY[0] = "Hemiechinus"
   }{
      print $1 "\t" LINARRAY[$3 + 0]
   }' > EH.popmap.txt
   rm admixture
fi

if [ ! -e complete.freq ]; then
   gawk --bignum \
        -v POPS="E_concolor,E_roumanicus_Balkan,E_roumanicus_Asian,E_europaeus_Iberian,E_europaeus_Apennine,Hemiechinus" \
        -v MINS="1,1,1,1,1,1" \
        -f $FREQ3 EH.popmap.txt EH.vcf > complete.freq
fi

# The file 'complete.freq', created with freq3.awk, contains the presumably
# derived allele's frequency in each of the given populations. The freq3.awk
# script does not assume a typical 4-taxon tree, but accepts any number of
# populations and calculates the derived allele frequencies in all them, in
# every site with at least a customizable number of individuals genotyped.
# One in this case. I have excluded E. amurensis in this file in order to
# optimize the number of sites with data. There are 261316 sites with at least
# one genotyped sample per population. Had I requested E. amurensis to also
# have one genotyped sample, I would have found less than 226000 sites.
#
# Among the six populations, how many site configurations are possible? Taking
# into account that the outgroup (Hemiechinus) is required to be fixed for the
# presumably ancestral allele, and distinguishing only ancestral and derived,
# there are 2^5=32 possible configurations: AAAAA, AAAAB, AAABA... BBBBB. Having
# over 260000 sites, it is feasible to estimate all their frequencies.
