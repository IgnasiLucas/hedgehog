2016-11-29
==========

It may be a good idea to take a look at the original bam files and use
bedtools to get an idea of the number of loci covered per individual and
the overlap among them.


2016-11-24
==========

The same conclusions than before are confirmed when including triallelic
variants in some statistics. The big picture is one of large genetic distance
between populations, and little variation within a population. This seems
to be the reason why little genetic structure is observed.

The question arises if the little genetic variation observed within
populations is due to a waste of reads in loci that don't get covered in
all samples, or if loci are quite homogeneously covered but turn out to be
actually invariant in many cases.

To answer this question, both variant and non-variant sites should have been
outputed by freebayes. In any case, the decay of the number of sites covered
in X samples when X increases suggests that the problem is one of too many
loci for not so many reads.


2016-11-18
==========

Using freebayes to call variants, and vcftools to summarize them, we
notice:
   1. Several types of variants, including SNP, indels, MNP, complex
      variants, and multiple combinations thereof when variable sites
      include more than two alleles.
   2. The vast majority of variable sites are fixed variants among
      populations, rather than variation segregating within populations.
