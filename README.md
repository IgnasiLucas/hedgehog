2018-07-25
==========

On 2018-07-24, Kristyna repeated the filtering of the original vcf file.
The new filtered vcf file should be the common starting point of all
downstream analysis. Here, I just add the BPF label to the INFO field
of that vcf file.

2018-07-19
==========

When trying to compare the SNPs identified as having europaeus ancestry
in the hybrid individual with the sites that show signal of historic
introgression between E. roumanicus and E. europaeus, Kristyna realized
that very few of those SNPs were included in the abba/baba analysis. In
this folder, I find out that most of those sites were arbitrarily excluded
from the abba/baba test because they were too close to other SNPs.


2018-06-30
==========

I wrote the script AncestryProfile.py, which takes the output files from
Admixture and uses an empirical Bayes approach to estimate the posterior
probability of each locus to have at least one allele coming from each
ancestry group. I run this script on all individuals, although I am only
interested in the hybrid between E. europaeus and E. romanicus. The script
seems to work well. Most of the genome show levels of ancestry probability
equal to the prior, which is the genome-wide contribution estimated by
Admixture. But highly informative SNPs produce spikes in the profile.

The main result is that the ~8% of E. europaeus ancestry in the hybrid
individual, is not clustered in discrete genomic blocks. Rather, the SNPs
with evidence of europaeus ancestry are scattered along the genome, very
mixed with romanicus-specific alleles. This suggests extensive recombination
and a more complex family history of this hybrid individual than just a
third generation of back-crossing.

2018-04-13
==========

I run ipyrad with the whole dataset. The results are available, but not
further used, yet.

2018-03-27b
===========

Preparation of a vcf file for abba/baba test. I want to optimize
the selection of sites so that there is enough data to perform the
abba/baba test. The strategy is not to rely on the general filters
available in vcftools, but to select sites with data for: the Hemichinus
sample (required for the test), and for any number of individuals from
the other three populations.

The script freq2.awk takes a vcf file with the "binary presence flag"
field in the INFO section (see 2017-05-25) as input and produces
a table of allele frequencies in the requested populations for
abba/baba analysis.

The filtering of the vcf file performed in this folder is superseded
by a common filtering strategy applied to the original vcf for all
downstream analyses. Thus, the thin.recode.vcf and all .tsv results
were erased from this folder. See 2018-07-25.

2018-03-27
==========

We received more sequence data, to improve the coverage. I should
re-run the alignments and ipyrad analysis. Kristýna has already
updated a vcf file. Here I just want to compare the new vcf file
to the old one, to evaluate how much more complete the data matrix
is after sequencing more.

2017-06-09
==========

Here I want to prepare the analysis with the multispecies coalescent
model, implemented in RevBayes. I need to format the data from the ipyrad
analysis and find out the most appropriate base substitution models for
each locus.

This analysis is interrupted because apparently RevBayes does not support
a too long list of loci. It is a good idea to filter loci, anyways,
because a Bayesian analysis will last forever if I use too many of them.

2017-05-25
==========

The current bam files were created by Kristýna with bowtie. They miss
read group information, but they are sorted by individual. I should add
the RG labels and run the freebayes pipeline to compare the results from
this bam file with those obtained before with a bowtie2-based alignment
(2016-11-18).

The results confirm that the data matrix is quite incomplete, to make any
analysis based on allele frequencies difficult. In this folder I add a
binary flag to the information field of the vcf file to indicate in what
samples each allele is present. This makes it easy to filter the vcf file
by sites covered, for example, in at least 5 individuals from E. romanicus
and 5 from E. europaeus.

2017-05-18
==========

Running ipyrad on the demultiplexed fastq files, using the reference and
the de novo assembly of reads that do not map to the reference. ipyrad
uses the bwa aligner. Apparently, it mapped 129,659,615 reads (40.5%) and
used 8,667,728 (2.7%) further reads with the de novo assembly pipeline.

After requesting the fixing of a bug, the output included the vcf file.

These results were not used, and superseded, anyways, by 2018-04-13. I
erase them.

2017-05-17
==========

Counting total numbers of raw and mapped reads, per sample. Overall the two
sequencing runs produced 319,934,047 reads (on average 6,398,680 reads per
sample, with a range of 717,180 -- 14,805,528). In all, 163290876 (51%) reads
were mapped with bowtie. 

2017-05-09
==========

Testing different mapping options to determine if we can improve the
mapping success and the coverage. I confirm that bowtie2 is more sensitive,
and that higher sensitivity leads to a larger proportion of multiple hits
(ambiguous mapping), and a decrease in unmapped reads, to a lesser extent.
Trimming the low quality tail of reads decreases the uniqueness of their
mappings with bowtie2. Unfortunately, bowtie1 does not report the numbers
of multiply mapped reads with sensitive settings. The best compromise
between a large portion of unique mappings and a low proportion of unmapped
reads is achieved using bowtie2 with whole reads and 'conservative' settings:
47.22% unique mappings, 18.46% unmapped.

To check the possible reasons for the unmapped reads, I manually did the
following (neither the commands nor the results are saved):

* Extract fastq files for mapped and unmapped reads from one bam file (Er71).
* Run fastqc on both to compare the qualities.
* Filter unmapped fastq records by maximum number of expected errors (0.5).
* BLAST a random subsample of 50 unmapped, high-quality reads to the nr database.

I conclude that the quality of unmapped reads is just slighly lower towards
the end of the reads. About 50% of the unmapped reads actually hit one or more
hedgehog sequences (either Erinaceus europeus, as expected, or Atelerix
albiventris, probably due to database limitations). Most such reads find
several hits, suggesting that they were unmapped (with bowtie1) due to
ambiguous mapping. Only one of 50 sequences mapped to a totally different
organism: a nematode (a parasite?). Contamination is, at most, rare.

2016-11-30
==========

Using bedtools, I compare the bam files among them, and find the number
of sites covered by at least 3 reads in N samples. I see, for example that
157661 sites are covered with at least 3 reads in at least 40 samples, which
seems to be plenty of information. Let's say that at least 60 bases of every
read are good enough to detect a variant. If heterozygosity was 0.001, I
would expect less than 10000 variants per individual. I suspect Freebayes
is not doing a very good job at identifying variants.

2016-11-24
==========

The same conclusions than before are confirmed when including combined
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
   1. Several types of variants, including SNP, MNP, complex variants,
      and multiple combinations thereof. The current freebayes settings
      seem to have excluded indels.
   2. The vast majority of variable sites are fixed variants among
      populations, rather than variation segregating within populations.
