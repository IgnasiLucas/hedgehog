#!/bin/bash
#
#				2020-09-28
#				==========
#
# Here I want to run some estimates of the f statistics, from the ABBA/BABA
# analysis, in a systematic way and using all data available. The goal is to
# compare different ways of computing the statistic in order to distinguish
# the ancient from the recent admixture signal.
#
# Recall the species tree:
#
#   (((E. concolor, E. roumanicus), E. europaeus), Hemiechinus sp.)
#   (((     P1    ,      P2      ),      P3     ),    OUTGROUP    )
#           \            /               /               /
#            \          /      ?        /               /
#             \        /  <-------->   /               /
#              \      /               /               /
#               \____/               /               /
#                  \                /               /
#                   \              /               /
#                    \____________/               /
#                           \                    /
#                            \                  /
#                             \                /
#                              \______________/
#                                      |
#
# Now, recall the results from Admixture, when using all individuals from the
# two contact zones: both E. roumanicus and E. europaeus have two lineages. And
# the current contact zones involve only the Eastern lineage of E. europaeus
# but with alternative lineages of E. roumanicus in either contact zone: the
# Balkanic lineage in Central Europe and the Asian lineage in the Russian contact
# zone.
#
#                                        E. roumanicus
#                           _______________________________________
#
#                               Balkanic               Asian
#                           _______________________________________
#
#    E. europaeus, Western     No contact            No contact
#
#    E. europaeus, Eastern   Current contact      Current contact
#   _______________________________________________________________
#
#
# The f statistics attempt to measure the average proportion of an individual's
# genome that is admixed; that is, that traces its genealogy through the other
# species (I think). Using individual lineages of the species involved in admixture,
# I expect the differences to be informative of the relative importance of admixture
# in different moments. Specifically:
#
#  -------------------------------------------------------------------------------------------------
#     Estimate              P2                          P3            Potential origin of admixture
#  -------------------------------------------------------------------------------------------------
#   f(ErB, EeW)   E. roumanicus, Balkanic     E. europaeus, Western              Ancient
#   f(ErB, EeE)   E. roumanicus, Balkanic     E. europaeus, Eastern      Ancient + recent in C.E.
#   f(ErA, EeW)   E. roumanicus, Asian        E. europaeus, Western              Ancient
#   f(ErA, EeE)   E. roumanicus, Asian        E. europaeus, Eastern      Ancient + recent in Russia
#  -------------------------------------------------------------------------------------------------
#
# Thus a measure of only ancient (previous to current interglacial) admixture would
# be the f statistic estimated with only the Western lineage of E. europaeus and the
# whole E. roumanicus. But, getting individual estimates f(ErB, EeW) and f(ErA, EeW)
# could reveal significant differences. I would expect f(ErB, EeW) > f(ErA, EeW), just
# because ancestors of the Balkanic E. roumanicus may have had the chance to meet
# ancestors of Western E. europaeus in more recent, previous interglacials than
# ancestors of Asian E. roumanicus have.
#
# The correct way to use the information would probably be to look directly at the
# frequencies of site patterns, and maybe using invariants to design statistics to
# test our hypotheses (Kubatko & Chifman, 2019, BMC Evolutionary Biology 19:112).
# But for the moment, the estimation of the fraction of introgression in 4-taxa
# subsets should be enough. In order to assess the uncertainty, I want to use a
# jackknive approach. I wonder whether deleting one individual at a time, instead
# of one genomic portion, would represent noise better.
#
# =================================================================================
#  SET UP
# =================================================================================
#
# System-specific variables
# -------------------------
#
# Awk script to add binary presence flag to a vcf file:
ADDBPF='../../bin/add_flag.awk'
FREQ2='../../bin/freq2.awk'
ADMIXTURE_DIR='/data/kristyna/hedgehog/results_2020/2020-06-12'
POPMAP_DIR='/data/kristyna/hedgehog/results_2020/2020-04-20'
SOURCE_VCF='/data/kristyna/hedgehog/results_2020/2020-02-28/all_merged_2020.vcf'
DATADIR=$(pwd | sed 's/results/data/')
LINEAGE=( "E_roumanicus_Asian" "E_europaeus_Iberian" "E_concolor" \
          "E_amurensis" "E_roumanicus_Balkan" "E_europaeus_Apennine" )
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
# clade of the E. europaeus species, wihtout any expected hybridization with E. roumanicus.

if [ ! -e EH.vcf ]; then
   if [ ! -e EH_list.txt ]; then
      cat $DATADIR/admixture.names.txt > EH_list.txt
      echo "Er65_IS25" >> EH_list.txt
   fi
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

# I decided to count patterns (BABA, ABBA) using only one individual per population.
# That way we get several measures, and an idea of the dispersion, even though the
# different measures are not completely independent, because the same individuals
# will be used more than once. The technical problem with this approach is that I may
# have to read the VCF files many times, one for every selection of population-specific
# individuals. For a given assignment of populations P1, P2, P3 and outgroup, if there
# are N1, N2 and N3 individuals from P1, P2 and P3, respectively, there will be
# N1*N2*N3 different combinations of one individual per population. When using the
# 24 Asian E. roumanicus, the 19 Apennine E. europaeus, and the 5 E. concolor, there
# will be 2280 combinations.
#
# Before anything else, let's use the freq2.awk script to get the derived allele
# frequencies for D and f estimation. The trees of interest are:
#
# (((E_concolor, E_roumanicus_Balkan), E_europaeus_Iberian), Hemiechinus)
# (((E_concolor, E_roumanicus_Balkan), E_europaeus_Apennine), Hemiechinus)
# (((E_concolor, E_roumanicus_Balkan), E_amurensis), Hemiechinus)
# (((E_concolor, E_roumanicus_Asian), E_europaeus_Iberian), Hemiechinus)
# (((E_concolor, E_roumanicus_Asian), E_europaeus_Apennine), Hemiechinus)
# (((E_concolor, E_roumanicus_Asian), E_amurensis), Hemiechinus)

if [ ! -d balkan_iberian ]; then mkdir balkan_iberian; fi
if [ ! -e balkan_iberian/complete.freq ]; then
   gawk --bignum \
        -v P1="E_concolor" \
        -v P2="E_roumanicus_Balkan" \
        -v P3="E_europaeus_Iberian" \
        -v OUTGROUP="Hemiechinus" \
        -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MINOUT=1 \
        -f $FREQ2 EH.popmap.txt EH.vcf | grep -v -F 999.9999 > balkan_iberian/complete.freq &
fi

if [ ! -d balkan_apennine ]; then mkdir balkan_apennine; fi
if [ ! -e balkan_apennine/complete.freq ]; then
   gawk --bignum \
        -v P1="E_concolor" \
        -v P2="E_roumanicus_Balkan" \
        -v P3="E_europaeus_Apennine" \
        -v OUTGROUP="Hemiechinus" \
        -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MINOUT=1 \
        -f $FREQ2 EH.popmap.txt EH.vcf | grep -v -F 999.9999 > balkan_apennine/complete.freq &
fi

if [ ! -d balkan_amurensis ]; then mkdir balkan_amurensis; fi
if [ ! -e balkan_amurensis/complete.freq ]; then
   gawk --bignum \
        -v P1="E_concolor" \
        -v P2="E_roumanicus_Balkan" \
        -v P3="E_amurensis" \
        -v OUTGROUP="Hemiechinus" \
        -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MINOUT=1 \
        -f $FREQ2 EH.popmap.txt EH.vcf | grep -v -F 999.9999 > balkan_amurensis/complete.freq &
fi

if [ ! -d asian_iberian ]; then mkdir asian_iberian; fi
if [ ! -e asian_iberian/complete.freq ]; then
   gawk --bignum \
        -v P1="E_concolor" \
        -v P2="E_roumanicus_Asian" \
        -v P3="E_europaeus_Iberian" \
        -v OUTGROUP="Hemiechinus" \
        -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MINOUT=1 \
        -f $FREQ2 EH.popmap.txt EH.vcf | grep -v -F 999.9999 > asian_iberian/complete.freq &
fi

if [ ! -d asian_apennine ]; then mkdir asian_apennine; fi
if [ ! -e asian_apennine/complete.freq ]; then
   gawk --bignum \
        -v P1="E_concolor" \
        -v P2="E_roumanicus_Asian" \
        -v P3="E_europaeus_Apennine" \
        -v OUTGROUP="Hemiechinus" \
        -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MINOUT=1 \
        -f $FREQ2 EH.popmap.txt EH.vcf | grep -v -F 999.9999 > asian_apennine/complete.freq &
fi

if [ ! -d asian_amurensis ]; then mkdir asian_amurensis; fi
if [ ! -e asian_amurensis/complete.freq ]; then
   gawk --bignum \
        -v P1="E_concolor" \
        -v P2="E_roumanicus_Asian" \
        -v P3="E_amurensis" \
        -v OUTGROUP="Hemiechinus" \
        -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MINOUT=1 \
        -f $FREQ2 EH.popmap.txt EH.vcf | grep -v -F 999.9999 > asian_amurensis/complete.freq &
fi

wait

# Note the use of --bignum flag in the gawk commands above. It is absolutely
# necessary for the large number of indiviudals present in the vcf file.
#
# The freq2.awk script trusts the binary presence flags to determine: what
# allele is ancestral (fixed in the outgroup, if any is fixed), and what
# the derived allele frequencies are in each population, in every valid site.
# The original D and f statistics were defined for a sample of one individual
# per population, and they are computed from genome-wide counts of site patterns
# (ABBA and BABA). When we have more than one individual per population,
# it is common practice to:
#
#   "... use the allele frequencies at each site to quantify the extent to which
#    the genealogy is skewed toward the ABBA or BABA pattern. This is effectively
#    equivalent to counting ABBA and BABA SNPs using all possible sets of four
#    haploid samples at each site. ABBA and BABA are therefore no longer binary
#    states, but rather numbers between 0 and 1 that represent the frequency of
#    allele combinations matching each genealogy. They are computed based on the
#    frequency of the derived allele (p) and ancestral allele (1- p) in each
#    population as follows:
#
#     ABBA = (1- p1 ) x p2 x p3 x (1- pO )
#     BABA = p1 x (1- p2 ) x p3 x (1- pO )"
#
#                                                                 Simon Martin, 2018
#                    http://evomics.org/learning/population-and-speciation-genomics/
#                      2018-population-and-speciation-genomics/abba-baba-statistics/
#
# However, I see a problem with this approach. By computing the probability of
# the ABBA and BABA patterns at every site as if all combinations of haplotypes
# were equally likely, we are assuming linkage equilibrium among sites. Not only
# linkage equilibrium is an unnecessary assumption for ABBA/BABA statistics, but
# it is known to be unlikely, to the point that the significance of the D statistic
# is assessed by jackknife with large genomic blocks precisely to take potential
# linkage disequilibrium into account (Durand et al. 2011, Mol. Biol. Evol. 28(8):
# 2239-2252). In the presenta analysis we are including introgressed or hybrid
# individuals, which make linkage disequilibrium a matter of fact. I have the
# impression the block jackknife procedure must underestimate the dispersion of
# the D and f statistics.
#
# Before attempting to count site patterns in all combinations of one individual
# per population, I will assess the dispersion by re-sampling the individuals
# used in the test, within populations. The outgroup species has only one sample.
# Population 1, E. concolor, is represented by 5 individuals, and I use sites
# covered in at least 1 of them, to optimize the number of sites. I consider
# the 5 pseudoreplicates produced by removing (or censoring) one of them at a
# time.
#
# The Iberian lineage of E. europaeus also has only 5 individuals, and I should
# use the same approach. E. amurensis is represented by only 2 individuals,
# allowing for 2 pseudoreplicates. All other populations have at least 19
# individuals, and I think it is more conservative to delete at least 5
# individuals at a time. That makes 42504 possible sub-samples of Asian
# lineage of E. roumanicus (with 24 individuals in all). And I should combine
# every subsample of P2 with every possible subsample of P3 (up to 11628,
# when using the Apennine lineage of E. europaeus, with 19 individuals).
# That means 2,471,182,560 combinations of all possible subsamples of
# E. concolor, E. roumanicus Asian and E. europaeus Apennine, with such
# a sub-sampling schema. Too many...

NUMREP=100

for pop in concolor balkan asian iberian apennine amurensis; do
   if [ ! -d $pop ]; then mkdir $pop; fi
done

N=1
for i in $(grep concolor EH.popmap.txt | cut -f 1); do
   if [ ! -e concolor/rep$N.txt ]; then
      echo $i > concolor/rep$N.txt
   fi
   N=$(( N + 1 ))
done

N=1
for i in $(grep amurensis EH.popmap.txt | cut -f 1); do
   if [ ! -e amurensis/rep$N.txt ]; then
      echo $i > amurensis/rep$N.txt
   fi
   N=$(( N + 1 ))
done

N=1
for i in $(grep Iberian EH.popmap.txt | cut -f 1); do
   if [ ! -e iberian/rep$N.txt ]; then
      echo $i > iberian/rep$N.txt
   fi
   N=$(( N + 1 ))
done

N=1
while [ $(md5sum apennine/* | cut -f 1 | sort | uniq | wc -l) -lt $NUMREP ]; do
   grep Apennine EH.popmap.txt | cut -f 1 | shuf -n 5 > $(printf "apennine/rep%03u.txt" $N)
   N=$(( N + 1 ))
done

N=1
while [ $(md5sum balkan/* | cut -f 1 | sort | uniq | wc -l) -lt $NUMREP ]; do
   grep Balkan EH.popmap.txt | cut -f 1 | shuf -n 5 > $(printf "balkan/rep%03u.txt" $N)
   N=$(( N + 1 ))
done

N=1
while [ $(md5sum asian/* | cut -f 1 | sort | uniq | wc -l) -lt $NUMREP ]; do
   grep Asian EH.popmap.txt | cut -f 1 | shuf -n 5 > $(printf "asian/rep%03u.txt" $N)
   N=$(( N + 1 ))
done

for P2 in balkan asian; do
   for P3 in iberian apennine amurensis; do
      OUTDIR=$(printf "%s_%s" $P2 $P3)
      N=1
      for con in $(ls -1 concolor); do
         for rou in $(ls -1 $P2); do
            for eur in $(ls -1 $P3); do
               cat concolor/$con $P2/$rou $P3/$eur > $OUTDIR/z1
               grep -v -F -f $OUTDIR/z1 EH.popmap.txt     > $(printf "%s/popmap%05u.txt" $OUTDIR $N)
               gawk '{print $1 "\tcensored"}' $OUTDIR/z1 >> $(printf "%s/popmap%05u.txt" $OUTDIR $N)
               N=$(( N + 1 ))
            done
         done
      done
      THREADS=1                                  # balkan_amurensis & asian_amurensis, 1000 replicates.
      if [ $N -gt 2000 ]; then THREADS=3; fi     # balkan_iberian & asian_iberian, 2500 replicates.
      if [ $N -gt 10000 ]; then THREADS=27; fi   # balkan_apennine & asian_apennine, 50000 replicates.
      PERTHREAD=$(( N / THREADS ))
      FIRST=1
      LAST=$PERTHREAD
      for i in $(seq 1 $THREADS); do
         printf "%s\t%s\t%u\t%u\t%u\n" $P2 $P3 $FIRST $LAST $N
         ./replicates.sh $P2 $P3 $FIRST $LAST &> $(printf "%s/rep%u_%u.log" $OUTDIR $FIRST $LAST) &
         FIRST=$(( FIRST + PERTHREAD ))
         LAST=$(( LAST + PERTHREAD))
      done
   done
done

wait
