2018-12-15
==========
It was clear from previous analyses that introgressed variation is not distributed
randomly along the genome in E. roumanicus. In order to identify introgressed variants
that have a higher (or lower) chance of lingering, Barbora suggested to estimate
genomic clines, using the Bayesian approach implemented in bgc by Gompert and Buerkle
(2012). For that, I first identify six individuals with varying degrees of admixture
between E. roumanicus and E. europaeus. Then, I run bgc using those six individuals
as the admixed population, and the rest as parental populations of either species.
Unfortunately, the admixed population is too small to notice any deviation of introgressed
allele frequencies with respect to the expectation of the genome-wide hybrid index.

2018-10-29
==========
Here I use Simon H. Martin's scripts to calculate population genetics statistics in
sliding windows. I also use the annotation of the reference genome to estimate those
statistics separately in genic and intergenic regions. I create some interesting plots
and found an unexpected negative relationship between diversity in E. roumanicus and
divergence between E. roumanicus and E. europaeus. The relationship is stronger in genic
than in intergenic regions, suggesting a functional origin (positive selection?).

2018-10-24
==========
Here I repeat the analysis from 2018-10-18, but with a different admixed individual,
Er55_AU7, which has only a small percentage of E. europaeus ancestry. Thus, its specific
combination of admixed and non-admixed genomic regions has actually survived more
generations and is more likely to be the result of natural selection.

  | Statistic |    Admixed         |   Not_admixed      |   Difference  |  p-value |
  | --------- | ------------------ | ------------------ | ------------- | -------- |
  | D         | 0.372100543001792  | 0.133008162873339  |    0.239092   |  0.06178 |
  | f_hom     | 0.106746542696069  | 0.0337681212184294 |    0.0729784  |  0.02096 |
  | f_d       | 0.0593813125368232 | 0.0198639588674299 |    0.0395174  |  0.05498 |
  | f         | 0.128633212368533  | 0.0388844749086059 |    0.0897487  |  0.02396 |

2018-10-18
==========
I compare the D and f statistics between genomic regions where the hybrid individual,
Er37_SK27, either has or does not have mixed ancestries. Both D and f statistics are
slightly higher where Er37_SK27 is admixed, even if Er37_SK27 itself is not used in
the computation of the statistics. However, differences are not significant, as assessed
by 5000 random permutations of the genomic regions. Here is the summary:

  | Statistic |  Admixed regions   |    Not admixed     |   Difference  |  p-value |
  | --------- | ------------------ | ------------------ | ------------- | -------- |
  | D         | 0.159047740296386  | 0.133996097982671  |   0.0250516   |   0.2678 |
  | f_hom     | 0.0398418118213496 | 0.0349892634711391 |   0.00485255  |   0.3332 |
  | f_d       | 0.0235819111306138 | 0.020606506731098  |   0.002975    |   0.3432 |
  | f         | 0.0452036938062747 | 0.0410208465643631 |   0.00418285  |   0.3542 |

2018-10-15
==========
I run the abba/baba test and the estimation of the proportion of the genome that
is introgressed. I use as input the most recent filtering of the vcf file, including
data from at least one individual of each population. I also split the E. europaeus
spcies in two populations, according to Admixture (and guessing the origin of two
individuals included in this analysis that were absent in the Admixture one). I
think I mostly confirm Kristýna's results. The D statistic, measuring the excess
of the ABBA pattern over the BABA one, is positive (0.15) and highly significant,
meaning that there must have been introgression between E. roumanicus and E. europaeus.
The signal is still strong when using only either the eastern or the western population
of E. europaeus. The proportion of the genome that is introgressed is low (below 4%,
and probably below 1%), and the estimate depends on the specific statistic used.

2018-10-08
==========
Simple update of the R script that creates the graphs from the Admixture output.
The goal was simply to order individuals in a reasonable way.

2018-10-03
==========
In the past (2018-06-30) I have used my own script to identify genomic regions
of admixture in individuals known to be admixed. However, my script did not take
recombination into account, as if all loci were independent. There exist several
programs that take recombination into account. The oldest of them do not require
phased haplotypes. I was unable to install SABER. I use LAMP instead. Some contigs
did not produce results (segmentation fault), but most did. I run LAMP with several
values for the time-since-admixture parameter, because for 3 of the 4 putatively
admixed individuals we don't know when admixture happened. Results are quite
consistent for a large range of those values. LAMP runs separately for each contig
and produces an image with ancestries represented as colors along a chromosome.
I joined the images corresponding to the 50 largest contigs to visualize the
most important results.


2018-09-26
==========

To address the problem detected on 2018-08-01, namely, the low accuracy of
allele frequency estimates in E. concolor population, I re-do the filtering
of the vcf file. Using our own scripts, we require either 5 or 4 E. concolor
individuals to have data. This reduces the number of sites to 43301 or 135043,
respectively, when also requiring at least 10 genotypes from both E. roumanicus
and E. europaeus. Kristýna uses these filtered vcf files to run again the
Admixture analysis. She finds that the small introgression signal from E. concolor
into (now two) individuals of E. roumanicus did not disappear.

2018-08-09
==========

Here, I check the effect of PCR duplicates. I count the number of second
reads that map in the same place, in all libraries; I look at the allele
count balance in individual libraries, and check Hardy-Weinberg statistics
The conclusion is that PCR duplicates are actually removing heterozygous
genotypes from the populations. In all, 16% sites show excess of homozygosity
in at least one population. However, allele frequencies should not be biased.

2018-08-01
==========

With the new filtering strategy, we find that one E. roumanicus individual
shows evidence of a small contribution by E. concolor. Introgression from
E. concolor is extremely unlikely, and the signal may be due to residual
incomplete lineage sorting (Barbora). Here, I look at what sites are giving
that signal, and make sure it is not an artifact. I confirm that there are
a few thousand sites, scattered across contigs where the presumably admixed
individual is heterozygous while the rest of E. roumanicus and E. concolor
specimens are (almost) fixed for alternative alleles. I also notice that many
of those sites do not include data for all 5 E. concolor individuals, and
I hypthesize that the lack of accuracy in allele frequency estimation there
is contributing to the signal.

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
