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
# allele. That is the data. Equation 2 in Alexander et al. (2009) gives the
# loglikelihood function. The likelihood of genotype g at one locus and
# one individual is:
#
#                         2!
# L(Q,P) = P(g|Q,P) = ---------- · [SUM q_k · p_k]^g · [SUM q_k · (1 - p_k)]^(2-g)
#                      g!(2-g)!      k                   k
#
# That is, the probability of finding g major alleles in a genotype equals
# the binomial probability of g successes out of 2, when the probability
# of an individual success is the probability of inheriting one such allele
# from any ancestry.
#
# I am interested in the probability of a specific ancestry
# in each locus, given the data and the major allele frequency in each ancestry,
# P.
