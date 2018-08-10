#!/bin/bash
#
#				2018-06-30
#				==========
#
# On 2018-06-05, Kristyna run Admixture to estimate the contributions of
# the different ancestral populations to each individual, as well as the
# allele frequencies in the ancestral populations. The dataset consisted
# on 67459 SNPs from a filtered vcf file, and 40 individuals. Howevere,
# there were some missing values. In one of the input files to the admixture
# program, I counted 2252611 genotypes available and 445749 missing genotypes.
# The number of parameters estimated is 67459 * K allele frequencies and
# K * 40 contributions from each ancestral population to each individual.
# For K = 4, which I think was the most appropriate one, this means there
# were 269996 parameters.
#
# Later on, we changed the filtering strategy of the vcf file, being more
# stringent in the quality of the variable site and the minimum distance
# between them, but more permissive in the quality of individual genotypes.
# This prevented many genotypes from being ignored. See Kristyna's folder
# 24-07-2018 for details. Currently, as 2018-07-26, the Admixture analysis
# in Kristyna's folder 05-06-2018 was run with 316489 variable sites, and
# 41 individuals. That makes almost 13 million genotypes (12976049), of
# which 15% (1942947) are missing.
#
# Here, I want to estimate the posterior probability of each ancestry in
# each locus, using as priors the genome-wide contributions of each ancestry
# to each individual. I think this is a naïve empirical Bayes approach,
# which I do not expect to be very accurate. I know I am ignoring the error
# in the estimates that I use as priors. The idea is simple: once we know
# the hybrid individual has a certain portion of E. europaeus ancestry,
# I would like to know what are the loci in its genome that more probably
# come from that ancestry.
#
# I am interested in the probability of each ancestry at each locus. Note
# that chromosomes may have different ancestries. Thus, I could estimate
# either the probability of ancestry k contributing to either copy of the
# locus, or the expected number of copies (between 0 and 2) contributed
# by the ancestry.
#
# Let's divide genotype g in homologous copies g0 and g1, each one being
# either the major or the minor allele. The probability of g0 coming from
# ancestry k is:
#
#               P(g0 | k0) · P(k0)
# P(k0 | g0) = --------------------
#                     P(g0)
#
# Where P(g0 | k0) is the frequency of allele g0 in ancestral population k,
# P(k0) is just the estimated contribution of k to the current individual,
# and P(g0) is the total probability of allele g0, across the K ancestries.
#
# The expected number of locus copies contributed by ancestry k would be
# P(k0 | g0) + P(k1 | g1). Otherwise, the overall probability that ancestry
# k contributed at least one allele to the genotype would be:
#
# P(k0 | g0) · P(not k1 | g1) + P(not k0 | g0) · P(k1 | g1) + P(k0 | g0) · P(k1 | g1)
#
# All of these can be calculated once I calculate the P(k0 | g0) and P(k1 | g1)
# for each k. For K=4, that is 8 probabilities per locus. For this I need
# The genotypes, the P and the Q matrices. While this could be done in R, I
# believe python will be better. The result should be a file like the P matrix,
# with as many lines as loci, with 2K columns for the alelle-specific ancestry
# probabilities, K columns for the expected number of alleles contributed by
# each ancestry, K columns for the probability of at least one allele contributed
# by each ancestry, and any additional columns. For example, the log odds
# between each ancestry and the next best one would make K additional columns
# useful for graphical representations. I will include contig and position.
# And all these, for each individual.
#
# When the observed allele frequency is 0.5, I don't know which allele's frequency
# Admixture reports. This is like not knowing what the genotype is, out of two
# options. I could calculate the two sets of posterior ancestry probabilities
# corresponding to the two possible genotypes, and then choose the genotype that
# makes the most likely overall ancestry of that individual (according to Admixture)
# also the most probable one a posteriori in that particular site.

ADMIXTURE_DIR=/data/kristyna/hedgehog/results_2018/05-06-2018
POPMAP_FILE=/data/joiglu/hedgehog/results/2018-03-27b/popmap.txt

if [ ! -e genotypes_vcfOrder.012 ]; then
   ln -s $ADMIXTURE_DIR/out.012 genotypes_vcfOrder.012
fi

if [ ! -e K3.P ]; then
   ln -s $ADMIXTURE_DIR/erinaceus_41.3.P K3.P
fi

if [ ! -e K3.Q ]; then
   ln -s $ADMIXTURE_DIR/erinaceus_41.3.Q K3.Q
fi

if [ ! -e SampleNames.txt ]; then
   cut -f 1 $ADMIXTURE_DIR/erinaceus_41.ped > SampleNames.txt
fi

if [ ! -e positions.txt ]; then
   # Note that here it is essential to use the .bim file, and not the .map file, since
   # they have different orders, and the .bim one is the one specifying the positions of
   # the variable sites in the binary .bed file used by admixture.
   cut -f 2 $ADMIXTURE_DIR/erinaceus_41.bim | sed 's/:/\t/' > positions.txt
fi

# Apparently, plink messes up the order of the genomic positions, because it does not
# consider more than 26 chromosomes. Thus, the genotypes in the 012 output of vcftools
# are not ordered in the same way as those actually used by admixture. And the genotypes
# used by admixture are in a binary file, difficult to parse. The best option is to change
# the order of the columns of the .012 file to match those in the binary .bed file. Below,
# I add line numbers to the list of positions in the wrong order. Then, I sequentially
# search for the positions in the right order, and use their line number to append the
# corresponding genotypes to a new genotypes file. It's very slow!

if [ ! -e genotypes.012 ]; then
   if [ ! -e positions_vcfOrder.txt ]; then
      # I add a tab at the end of the line to match it later, when searching.
      nl $ADMIXTURE_DIR/out.012.pos | gawk '{print $1 "\t" $2 "\t" $3 "\t"}' > positions_vcfOrder.txt
   fi
   NUMLOCI=`cat positions.txt | wc -l`
   cut -f 1 genotypes_vcfOrder.012 > genotypes.012
   for locus in `seq 1 $NUMLOCI`; do
      head -n $locus positions.txt | tail -n 1 | gawk '{print $1 "\t" $2 "\t"}' > pattern.txt
      VCFORDER=`grep -Ff pattern.txt positions_vcfOrder.txt | cut -f 1`
      cut -f $(( $VCFORDER + 1 )) genotypes_vcfOrder.012 > z1
      paste genotypes.012 z1 > z2
      mv z2 genotypes.012
      rm z1 pattern.txt
   done
fi

if [ ! -e popmap.txt ]; then
   ln $POPMAP_FILE popmap.txt
fi

# Here, I want to see if the ancestry groups identified correspond to the
# known populations.
if [ ! -e summary_Q.txt ]; then
   echo "Expected number of individuals from each ancestry:" > summary_Q.txt
   gawk '{for (i = 1; i <= NF; i++) {N[i] += $i}}END{for (i in N) print i "\t" N[i]}' K3.Q | sort -nk 1,1 >> summary_Q.txt
   for k in 1 2 3; do
      for sampleIndex in `gawk -v K=$k '($K > 0.5){print NR}' K3.Q`; do
         head -n $sampleIndex SampleNames.txt | tail -n 1 >> z$k.txt
      done
      echo "Individuals with more than 50% ancestry $k:" >> summary_Q.txt
      grep -f z$k.txt popmap.txt >> summary_Q.txt
      rm z$k.txt
   done
fi
# With K=4, the program Admixture correctly identified E. concolor and E. romanicus.
# The species E. europaeus was split in two populations, a western and an eastern one,
# I think. From the expected number of individuals from each ancestry, the expected
# current allele frequency can be estimated for each locus in K4.P. I manually checked
# that the predicted allele frequencies are usually above 0.5. That is, Admixture is
# actually using the major allele to report allele frequencies, just as claimed in the
# paper. It is important to take into account that the file 'genotypes.012' is not
# expressed as number of major alleles, but probably as number of alternative alleles.
# That is, I need to translate the genotypes.

#for ind in `cat SampleNames.txt`; do
for ind in Er37_SK27 Er55_AU7 Er50_R3 Er26_JUG4; do
   if [ ! -e $ind.out ]; then
      python AncestryProfile.py -i $ind \
                                -s SampleNames.txt \
                                -k 3 \
                                -p K3.P \
                                -q K3.Q \
                                -g genotypes.012 \
                                -m positions.txt \
                                -o $ind.out
   fi
done

# I expect the alternative ancestries of a hybrid individual to be well
# separated in alternative blocks along the contigs. However, it does not
# need to be the case. Some SNPs should be very informative of the ancestry
# while others will not.
#
# There are 1747 contigs with data. Only 189 of them have at least 100 SNPs. I will
# plot only those for the hybrid individual. I get plots as well from one individual
# representative of each ancestry group, for comparison. I choose the individual
# in each group with the least number of missing genotype values.

MINIMUM=200   # Minimum number of SNPs in a contig to plot the ancestry profile.
for sample in Er37_SK27 Er55_AU7 Er50_R3 Er26_JUG4; do
   if [ ! -d $sample ]; then mkdir $sample; fi
   if [ ! -d $sample/LOD ]; then mkdir $sample/LOD; fi
   if [ ! -d $sample/ENum ]; then mkdir $sample/ENum; fi
   touch $sample/NotPlotted.txt
   for contig in `cut -f 1 $sample.out | sort | uniq`; do
      if [ ! -e $sample/$contig.png ] && ! grep -q $contig $sample/NotPlotted.txt; then
         NUM_SNPS=`grep $contig $sample.out | tee $sample/zdata.txt | wc -l`
         if [ $NUM_SNPS -ge $MINIMUM ]; then
            if [ ! -e $sample/LOD/$contig.png ]; then
               gnuplot -e "contig='$contig'; infile='$sample/zdata.txt'; outfile='$sample/LOD/$contig.png'" PlotAncestryProfile.gnp
            fi
            if [ ! -e $sample/ENum/$contig.png ]; then
               gnuplot -e "contig='$contig'; infile='$sample/zdata.txt'; outfile='$sample/ENum/$contig.png'" PlotExpectedAlleleNumber.gnp
            fi
         else
            if ! grep -q $contig $sample/NotPlotted.txt; then
               echo "Adding contig $contig to $sample/NotPlotted.txt."
               cat $sample/zdata.txt >> $sample/NotPlotted.txt
            fi
         fi
         rm $sample/zdata.txt
      fi
   done
done

# Finally, I want a list of loci where the contribution of ancestry group 1
# to the hybrid individual (with group 2 background under K=3) is probable enough.
# Say, with logodds > 3.0

MIN_LOGODDS=3
if [ ! -e Er37_SK27.europaeus_SNPs.txt ]; then
   gawk -v MINLOD=$MIN_LOGODDS '($3 >= MINLOD){print $0}' Er37_SK27.out > Er37_SK27.europaeus_SNPs.txt
fi


# CONCLUSIONS
# ===========
#
# The contribution of an europaeus ancestry to the genome of the hybrid
# individual is scattered along more than 2000 loci in 788 contigs. After
# inspecting visually the ancestry profiles of 189 contigs with at least
# 100 SNPs each, I have not seen a block of europaeus ancestry spanning
# more than two SNPs. About 81% of loci with europaeus ancestry also show
# signs of romanicus ancestry. This can be caused by either heterozygosity
# or by a similar estimated allele frequency in both ancestries.
#

if [ ! -e summary.txt ]; then
   echo -e "#Number of loci probably contributed by each ancestry (log odds > 3.0)" > summary.txt
   echo -e "#Sample\tAnc.1\tAnc.2\tAnc.3" >> summary.txt
   for sample in `cat SampleNames.txt`; do
      if [ -e $sample.out ]; then
         gawk -v SAMPLE=$sample '{
            for (i=3; i<=NF; i++) {
               if ($i > 3) N[i-2]++
            }
         }END{
            printf("%10s\t%i\t%i\t%i\t%i\n", SAMPLE, N[1] + 0, N[2] + 0, N[3] + 0, N[4] + 0)
         }' $sample.out >> z1
      fi
   done
   sort -nrk 4 -rk 2 -rk 5 -rk 3 z1 >> summary.txt
   rm z1
fi

# Overall, I understand this to mean that there are at least some alleles quite
# specific and informative of the specific ancestries.
