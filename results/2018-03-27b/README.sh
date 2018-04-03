#!/bin/bash
#
#				2018-03-27b
#				===========
#
# Kristýna has previously used a set of scripts to apply the abba/baba test,
# that are available here: https://github.com/simonhmartin/genomics_general
#
# The usage of those scripts is documented in the following website:
#    http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/
#
# The freq.py script parses a matrix of genotypes and outputs the derived
# allele frequencies for a set of sites, in each of three populations. The
# fourth population is the outgroup, which must be fixed for the ancestral
# population.
#
# The use of population-specific allele frequencies to run the abba/baba test
# is convenient for an optimal use of the data, because it allows us to use
# any (reasonable) number of individuals of each population in every site.
# That is, we do not need a complete matrix with the exact same individuals
# genotyped in all loci.
#
# However, the freq.py script is not optimal. In order to calculate allele
# frequencies from genotypes, it splits the genotypes in as many alleles as
# the ploidy number, and then counts the number of times each base is observed
# among all the observed chromosomes. It requires a relatively high confidence
# in the genotype, say at least 5 reads supporting it. Otherwise, the two
# chromosomes in a diploid may be undersampled. If we did not really sample
# 2n chromosomes, but somewhat less, the standard error of the frequency must
# be higher than expected. At first I also thought it could introduce some bias
# under some circumstances.
#
# In folder 'estimation', I run a simulation of reads drawn from a set of
# individuals with known allele and genotype frequencies. Then, I estimate the
# allele frequencies either assuming we know the genotypes with accuracy or
# applying a correction. As it can be easily shown analytically, the simplification
# of assuming that we know the genotypes with accuracy does not bias the allele
# frequency estimates. About the variance of the estimates, I did not do the math,
# but from the simulations it seems that it can be improved just slightly, for
# low coverage values.
#
# What I did not check is the effect of excluding individuals with low coverage.
#
# For the moment, I will focus on getting a vcf file and the associated frequency
# file.

GENOMICS_GENERAL=~/bin/genomics_general

if [ ! -d estimation ]; then mkdir estimation; fi
if [ ! -e estimation/variances.png ]; then
   R --save < simulate.R 1> estimation/log 2> estimation/err &
fi

VCF=/data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz

if [ ! -e popmap.txt ]; then
   cp ../2017-05-25/populations.txt ./popmap.txt
fi

if [ ! -e thinned.recode.vcf ]; then
   if [ ! -e filtered.vcf ]; then
      if [ ! -e flagged.vcf ]; then
         gunzip -c $VCF > z1.vcf
         gawk -f add_flag.awk z1.vcf > flagged.vcf
      fi
      rm z1.vcf
      gawk -v MINQ=50 -f filtervcf.awk popmap.txt flagged.vcf > filtered.vcf
   fi
   rm flagged.vcf
   vcftools --vcf filtered.vcf \
            --out thinned \
            --thin 200 \
            --recode \
            --recode-INFO NS \
            --recode-INFO AN \
            --recode-INFO AF \
            --recode-INFO ABP \
            --recode-INFO BPF
fi
rm filtered.vcf

# Now, before writing my own script to generate the derived allele frequencies,
# I will try to apply Simon Martin's pipeline to the thinned.recode.vcf file.
# I have downloaded the scripts in my home's bin directory. Change the path to
# genomics_general above to run somewhere else.

if [ ! -e all.geno.gz ]; then
   python $GENOMICS_GENERAL/VCF_processing/parseVCF.py -i thinned.recode.vcf \
                                                       --ploidy 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 \
                                                       --skipIndels | gzip > all.geno.gz
fi

if [ ! -e all1.tsv ]; then
   python $GENOMICS_GENERAL/freq.py -g all.geno.gz \
       -p rom Er26_JUG4,Er27_SK32,Er28_112,Er32_183,Er35_209,Er36_SK24,Er39_PL1,Er40_M2,Er41_GR36,Er42_GR35,Er43_SL7,Er44_VOJ1,Er45_BLG3,Er46_RMN7,Er47_CR4,Er48_BH16,Er50_R3,Er55_AU7,Er60_GR5,Er61_GR87,Er62_GR95,Er66_IT3,Er69_R2,Er70_RMN42 \
       -p eur Er29_122,Er30_79,Er31_453,Er33_211,Er34_197,Er51_436,Er52_451,Er54_AU1,Er56_AZ5,Er57_COR4,Er58_FI7,Er59_FR1,Er63_IR6,Er67_IT5,Er68_PRT1B,Er71_SAR2,Er72_SIE1B,Er74_SP16 \
       -p con Er38_LB1,Er49_GR38,Er53_ASR7,Er64_IS1,Er75_TRC2A \
       -p hem Er65_IS25 \
       --target derived \
       -o all1.tsv
fi

# I think my script is quite faster, although it can use only one thread. The MINX
# parameters below set the minimum number of individuals required per population.
# Populations 1 to 3 are E. romanicus, E. europaeus, and E. concolor. The outgroup
# is Hemiechinus. This is defined in the script itself.

if [ ! -e all2.tsv ]; then
   gawk -v MIN1=5 -v MIN2=5 -v MIN3=4 -v MINOUT=1 -f freq.awk popmap.txt thinned.recode.vcf > all2.tsv
fi

# I do not reproduce here a manual check: transforming the tsv files to bed files,
# and using bedtools to find the common positions, I checked that the frequencies
# calculated with my script are exactly the same as those calculated with Simon
# Martin's script. The sites reported are different, though, because of different
# filters.
#
# The main point from this analysis is that using individuals with low coverage,
# where undersampling of heterozygous genotypes may be likely, does not bias the
# estimate of allele frequency. Thus, as long as we trust the variable site, because
# it is variable in several populations, even individuals with only 1 read contribute
# valuable information about the allele frequency, and should not be excluded from
# the estimate. Even if the standard error of the estimated allele frequency is higher
# than it would be under high coverage, adding one read from one individual actually
# lowers the standard error, with respect to the exclusion of that individual.
