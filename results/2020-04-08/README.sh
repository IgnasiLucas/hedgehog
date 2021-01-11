#!/bin/bash
#
#				2020-04-08
#				==========
#
# The unexpected results in 2020-04-07 make me suspect that there could be
# an artifact. We are expecting a positive correlation between within-species
# diversity and between-species divergence, on the bases of well established
# neutral theory. I need to check how diversity and divergence is measured.
#
# The data is a collection of genotypes at SNP sites: loci selected because
# they are variable. Genotypes are turned into an alignment, with alleles
# of one diploid individual distributed in two sequences. The distance between
# two individuals is the average Hamming distance between their sequences,
# without counting positions with Ns. I don't see where in Martin's code, but
# at some point this distance should be normalized by the length of the alignment
# or number of sites used. Thus, nucleotide diversity, pi, in a species is
# the expected proportion of nucleotides that differ between two randomly
# chosen sequences from the same population, among the number of sites used.
# And the divergence, dxy, is the expected proportion of differences between
# sequences from two different populations.
#
# Could the ascertainment bias cause artifacts in the correlation between
# pi and dxy? We are using a SNP matrix with at least 10 valid genotypes in
# E. europaeus, 10 more in E. roumanicus, and 4 (out of 5) in E. concolor.
# We know from the tree of species that most sites show differences
# between E. europaeus and the other two. And the fact is that our power to
# detect variation is higher in E. roumanicus and E. europaeus than in
# E. concolor.
#
# The correlations appear when considering more than one locus, because of
# the variation in the distribution of the different kinds of variable sites
# along the genome. Thus, it appears that if a window contains many divergent
# sites between E. roumanicus and E. europaeus, then it usually contains relatively
# few sites variable within E. roumanicus.
#
# One possible artifactual reason for this is the fact of defining windows by
# the number of variable sites in them. If the number of variable sites in a
# window is limited, too many of one kind will necessarily meen fewer of the
# other kind. So, the proper way to define windows in this case is by coordinates,
# so that the window limits are random with respect to the distribution of the
# variable sites.
#
# The ascertainment bias could also be a problem. That is, the fact that sites
# that take part in the alignment are not a random sample of the genome: divergent
# sites and variable sites in different species have different probabilities of
# being included in the analysis.

# ---------------------------------------------------------------------------------
#  SET UP
# ---------------------------------------------------------------------------------

VCF=../2018-09-26/ErinMaxMiss63_r10e10c4.vcf
GENOMICS_GENERAL=~/bin/genomics_general
POPDATA=../../data/populations.txt
NUMREP=100

if [ ! -e erin63.geno.gz ]; then
   python $GENOMICS_GENERAL/VCF_processing/parseVCF.py -i $VCF | gzip > erin63.geno.gz
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

   if [ ! -e erin63.popmap.txt ]; then
      gunzip -c erin63.geno.gz | head -n 1 | sed 's/\t/\n/g' | tail -n +3 > good_samples.txt
      grep -F -f good_samples.txt popmap.txt > erin63.popmap.txt
      rm good_samples.txt
      rm popmap.txt
   fi
fi

# ---------------------------------------------------------------------------------
# WINDOW DEFFINITION
# ---------------------------------------------------------------------------------
#
# I will not care now about introgression statistics, in order to have more
# flexibility in the deffinition of window sizes. The goal in this section is
# to determine the effect of the genomic window size in the correlations between
# divergence and diversity.

if [ ! -d windowSize ]; then mkdir windowSize; fi
for W in 500000 750000 1000000 1250000 1500000 1750000 2000000 2250000 2500000; do
   if [ ! -e windowSize/$W.popGenStats.csv ]; then
      python $GENOMICS_GENERAL/popgenWindows.py --windType coordinate -w $W \
         --overlap 0 -m 30 -f phased -T 2 \
         --popsFile erin63.popmap.txt \
         -p concolor -p roumanicus -p europaeus \
         -g erin63.geno.gz \
         -o windowSize/$W.popGenStats.csv &
   fi
   if [ ! -e windowSize/$W.popFreq.csv ]; then
      python $GENOMICS_GENERAL/popgenWindows.py --windType coordinate -w $W \
         --overlap 0 -m 30 -f phased -T 2 \
         --popsFile erin63.popmap.txt \
         --analysis popFreq \
         -p concolor -p roumanicus -p europaeus \
         -g erin63.geno.gz \
         -o windowSize/$W.popFreq.csv &
   fi
done
wait

# --------------------------------------------------------------------------------
# ASCERTAINMENT BIAS
# --------------------------------------------------------------------------------
#
# Here, I want to subsample the populations with a larger number of individuals.
# That will make some previously variable sites invariable, therefore balancing
# the power to detect variation among populations.

if [ ! -d subsampling ]; then mkdir subsampling; fi

grep concolor   erin63.popmap.txt | cut -f 1 > concolor.txt
grep roumanicus erin63.popmap.txt | cut -f 1 > roumanicus.txt
grep europaeus  erin63.popmap.txt | cut -f 1 > europaeus.txt
gunzip -c erin63.geno.gz | head -n 1 | sed 's/\t/\n/g' | nl -n ln | grep -f concolor.txt   | cut -f 1 > concolor.columns.txt
gunzip -c erin63.geno.gz | head -n 1 | sed 's/\t/\n/g' | nl -n ln | grep -f roumanicus.txt | cut -f 1 > roumanicus.columns.txt
gunzip -c erin63.geno.gz | head -n 1 | sed 's/\t/\n/g' | nl -n ln | grep -f europaeus.txt  | cut -f 1 > europaeus.columns.txt
rm concolor.txt roumanicus.txt europaeus.txt

gawk '(NR == 1){S = $1}(NR > 1){S = S "," $1}END{print S}' concolor.columns.txt > concolor_combination.txt
touch roumanicus_combinations.txt
while [[ $(cat roumanicus_combinations.txt | wc -l) -lt $NUMREP ]]; do
   # Happy to find the shuf command!
   shuf roumanicus.columns.txt | head -n 5 | sort -n | gawk '(NR == 1){S = $1}(NR > 1){S = S "," $1}END{print S}' >> roumanicus_combinations.txt
   sort roumanicus_combinations.txt | uniq > z
   mv z roumanicus_combinations.txt
done
rm roumanicus.columns.txt

touch europaeus_combinations.txt
while [[ $(cat europaeus_combinations.txt | wc -l) -lt $NUMREP ]]; do
   shuf europaeus.columns.txt | head -n 5 | sort -n | gawk '(NR == 1){S = $1}(NR > 1){S = S "," $1}END{print S}' >> europaeus_combinations.txt
   sort europaeus_combinations.txt | uniq > z
   mv z europaeus_combinations.txt
done
rm europaeus.columns.txt

for i in $(seq 1 $NUMREP); do
   if [[ ! -e $(printf "subsampling/%03i.geno.gz" $i) ]]; then
      gunzip -c erin63.geno.gz | \
      cut -f 1,2,$(cat concolor_combination.txt),$(head -n $i roumanicus_combinations.txt | tail -n 1),$(head -n $i europaeus_combinations.txt | tail -n 1) | \
      gawk '(NR == 1){print}(NR > 1){
         delete P
         for (i = 3; i <= NF; i++) {
            split($i, Z, /\//)
            P[Z[1]] = 1
            P[Z[2]] = 1
         }
         delete P["N"]
         if (length(P) > 1) print
      }' | \
      gzip - > $(printf "subsampling/%03i.geno.gz" $i)
   fi
   if [ ! -e $(printf "subsampling/%03i.map.txt" $i) ]; then
      gunzip -c $(printf "subsampling/%03i.geno.gz" $i) | head -n 1 | sed 's/\t/\n/g' | tail -n +3 > z
      grep -f z erin63.popmap.txt > $(printf "subsampling/%03i.map.txt" $i)
      rm z
   fi
   if [ ! -e $(printf "subsampling/%03i.popGenStats.csv" $i) ]; then
      python $GENOMICS_GENERAL/popgenWindows.py --windType coordinate -w 750000 \
         --overlap 0 -m 30 -f phased -T 2 \
         --popsFile $(printf "subsampling/%03i.map.txt" $i) \
         -p concolor -p roumanicus -p europaeus \
         -g $(printf "subsampling/%03i.geno.gz" $i) \
         -o $(printf "subsampling/%03i.popGenStats.csv" $i)
   fi
   if [ ! -e $(printf "subsampling/%03i.popFreq.csv" $i) ]; then
      python $GENOMICS_GENERAL/popgenWindows.py --windType coordinate -w 750000 \
         --overlap 0 -m 30 -f phased -T 2 \
         --popsFile $(printf "subsampling/%03i.map.txt" $i) \
         --analysis popFreq \
         -p concolor -p roumanicus -p europaeus \
         -g $(printf "subsampling/%03i.geno.gz" $i) \
         -o $(printf "subsampling/%03i.popFreq.csv" $i)
   fi
done
rm concolor_combination.txt

if [ ! -e artifacts.pdf ]; then
   R -q --no-save -e "rmarkdown::render('artifacts.Rmd', output_file='artifacts.pdf')"
fi
