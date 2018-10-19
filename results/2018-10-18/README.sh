#!/bin/bash
#
#				2018-10-18
#				==========
#
# Once I am familiar with the D-test (2018-10-15) and I know what parts of the
# genome of the admixed Er37_SK27 individual have an E. europaeus origin (2018-10-03),
# it's time to determine if the genomic regions where Er37_SK27 keeps E. europaeus
# ancestry are also characterized by a strong introgression signal among its peers.
#
# The idea may be too simple. Er37_SK27 must be the offspring of one normal and
# one hybrid, fully heterozygous parent. Hence, its approximately 25% of E. europaeus
# ancestry. Because Er37_SK27 survived development, and to the moment of sample
# collection, it can be said that its specific combination of ancestyries is not
# lethal, at least. And, because it survived that long, the most likely scenario is
# that the combination of ancestries is not deleterious. If natural selection removes
# from the population deleterious admixed genomic regions, then we expect signatures
# of admixture to be localized only where natural selections allows it. If that was
# most of the genome, we would not expect much coincidence among independently admixed
# individuals. But, if admixture was deleterious in most of the genome, even individuals
# with low levels of admixture who are not related, would share the admixed blocks.
#
# So, what I want to test is if the D statistic in the admixed regions is significantly
# higher than that outside the admixture regions. To assess significance, I need a null
# distribution... But let's check the difference first.

if [ ! -e admixture_regions.tsv ]; then
   
fi

if [ ! -e not_admixed_regions.tsv ]; then

fi

if [ ! -e D.txt ]; then

fi
