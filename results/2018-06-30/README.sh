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
# program, I count 2252611 genotypes available and 445749 missing genotypes.
# The number of parameters estimated is 67459 * K allele frequencies and
# K * 40 contributions from each ancestral population to each individual.
# For K = 4, which I think was the most appropriate one, this means there
# are 269996 parameters.
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
# There is a file already available, with extension "012" and created by
# vcftools, where genotypes have been translated into counts of the major
# allele. That is the data.
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
#
# I took the following files from Kristyna's folder 05-06-2018, with minor
# modifications of names or format. K4.p and K4.Q are Admixture's output.
# genotypes.012 was named out.012. The SampleNames.txt is the first
# column of erinaceus.ped, and positions.txt is the second column of
# erinaceus.map (with ":" substituted by tab). The popmap.txt is from 2018-03-27b.

for file in K4.P K4.Q genotypes.012 SampleNames.txt positions.txt popmap.txt; do
   if [ ! -e $file ]; then
      echo "File $file not found."
      exit
   fi
done

# Here, I want to see if the ancestry groups identified correspond to the
# known populations.
if [ ! -e summary_Q.txt ]; then
   echo "Expected number of individuals from each ancestry:" > summary_Q.txt
   gawk '{for (i = 1; i <= NF; i++) {N[i] += $i}}END{for (i in N) print i "\t" N[i]}' K4.Q | sort -nk 1,1 >> summary_Q.txt
   for k in 1 2 3 4; do
      for sampleIndex in `gawk -v K=$k '($K > 0.5){print NR}' K4.Q`; do
         head -n $sampleIndex SampleNames.txt | tail -n 1 >> z$k.txt
      done
      echo "Individuals with more than 50% ancestry $k:" >> summary_Q.txt
      grep -f z$k.txt popmap.txt >> summary_Q.txt
      rm z$k.txt
   done
fi
# With K=4, the program Admixture correctly identified E. concolor and E. romanicus.
# The species E. europaeus is split in two pupulations, a western and an eastern one,
# I think. From the expected number of individuals from each ancestry, the expected
# current allele frequency can be estimated for each locus in K4.P. I manually checked
# that the predicted allele frequencies are usually above 0.5. That is, Admixture is
# actually using the major allele to report allele frequencies, just as claimed in the
# paper. It is important to take into account that the file 'genotypes.012' is not
# expressed as number of major alleles, but probably as number of alternative alleles.
# That is, I need to translate the genotypes.

for ind in `cat SampleNames.txt`; do
   if [ ! -e $ind.out ]; then
      python AncestryProfile.py -i $ind \
                                -s SampleNames.txt \
                                -k 4 \
                                -p K4.P \
                                -q K4.Q \
                                -g genotypes.012 \
                                -m positions.txt \
                                -o $ind.out
   fi
done

# I expect the alternative ancestries of a hybrid individual to be well
# separated in alternative blocks along the contigs. However, it does not
# need to be the case.
