#!/bin/bash
#
#				2017-06-09
#				----------
#
# Loci with more than two parsimony informative sites are easy to filter
# from the alleles.loci output of ipyrad. I should output each alignment
# to a different file. I think it is safe to skip the second allele of an
# individual if it is identical to the first one. I do not want too many
# sequences to analyse, and most such cases must indicat undersampling,
# rather than real homozygosity.
#
# I cannot use the nexus output file from ipyrad, because it concatenates
# all loci, assuming that all loci must be sampled from the same set of
# individuals. Actually, it is interesting to note that when keeping loci independent
# it does not matter if the chromosomes sampled at each locus are different.
#
# It is probably easy to parse the output alignments directly from within the
# API of ipyrad, using python commands to filter and output the alignments in
# different nexus files. But the API is not well documented, and I have no
# idea how to do that.

LOCI=../2017-05-18/min5_outfiles/min5.alleles.loci

if [ ! -d nexus ]; then mkdir nexus; fi

if ! ls nexus | grep -q ".nex"; then
   gawk -f loci2nex.awk $LOCI
   mv *.nex nexus/
fi
