#!/bin/bash
#
#				2018-10-29
#				==========
#
# After having shown that regions of the genome where Er55_AU7 retain E. europaeus
# ancestry from a relatively recent admixture event are also enriched in signatures
# of historical introgression among the rest of individuals (detected with the D
# test), it is necessary to characterize those regions. If introgressed variation
# lingers where it is less harmful, it would be found mostly where functional
# constraint is low. Functional constraint could be detected by looking at the
# annotation files. Also, comparing polymorphism and divergence between two or more
# loci (Hudson, Kreitman, Aguadé test), we can gain insight.
#
# Also, there is no need to limit our analysis to the regions where Er55_AU7 retained
# E. europaeus ancestry. That's just one case, and others would have improved our power.
# The question is that with that case we proved an association between two types of
# introgression signal. Now, we can look for an association between introgression
# signal and functional constraint. For example, estimating D and f statistics both
# for genic and intergenic regions. Also, estimating either D or f statistics in
# windows, where both divergence and polymorphism is also estimated. The latter can be
# done with Simon H. Martin's scripts: https://github.com/simonhmartin/genomics_general

FREQ2_AWK=../../bin/freq2.awk
# Different vcf files are optimized for different statistics:
VCF1=../2018-09-26/ErinaceusAndHemiechinus_r1e1c1h1.vcf
VCF2=../2018-09-26/ErinMaxMiss63_r10e10c4.vcf
# This is the folder where Simon H. Martin's genomic_general repository is cloned. It
# works with python 2, not 3. I add a conda env for this current folder. Everything is
# executed within that conda environment, with pyton 2, that I call '2018-10-29'
GENOMICS_GENERAL=~/bin/genomics_general
POPDATA=../../data/populations.txt

# I will exclude both the hybrid and the other admixed individual, in order to work only
# with the historical signal of introgression. It will help interpretation. The hybrid
# individual is already identified as such in the popmap file. Below, I assign Er55_AU7
# to a population called 'admixed'. That's enough for them to be ignored when running
# either popgenWindows.py or ABBABABAwindows.py.

if [ ! -e creh1.geno.gz ]; then
   python $GENOMICS_GENERAL/VCF_processing/parseVCF.py -i $VCF1 | gzip > creh1.geno.gz
fi


if [ ! -e erin63.geno.gz ]; then
   python $GENOMICS_GENERAL/VCF_processing/parseVCF.py -i $VCF2 | gzip > erin63.geno.gz
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

   if [ ! -e creh1.popmap.txt ]; then
      gunzip -c creh1.geno.gz | head -n 1 | sed 's/\t/\n/g' | tail -n +3 > good_samples.txt
      grep -F -f good_samples.txt popmap.txt > creh1.popmap.txt
      rm good_samples.txt
   fi
   if [ ! -e erin63.popmap.txt ]; then
      gunzip -c erin63.geno.gz | head -n 1 | sed 's/\t/\n/g' | tail -n +3 > good_samples.txt
      grep -F -f good_samples.txt popmap.txt > erin63.popmap.txt
      rm good_samples.txt
      rm popmap.txt
   fi
fi

# The abbababa statistics are quite demanding on the number of SNPs used, and they are necessarily
# based on a smaller portion of the original VCF file, because they require data for the one Hemiechinus
# sample. Thus, I'll run first the abbababa tests in sliding windows and use those windows coordinates
# to calculate population genetic statistics with the larger dataset.

if [ ! -e creh1.abbababa.csv ]; then
   python $GENOMICS_GENERAL/ABBABABAwindows.py -g creh1.geno.gz -f phased -o creh1.abbababa.csv --windType sites \
      -w 110 --overlap 100 -m 100 --popsFile creh1.popmap.txt -P1 concolor -P2 roumanicus -P3 europaeus -O Hemiechinus --writeFailedWindows -T 20
fi

if [ ! -e erin63.PopGenStats.csv ]; then
   if [ ! -e creh1.windows.txt ]; then
      cut -f 1,2,3 -d ',' --output-delimiter=$'\t' creh1.abbababa.csv | tail -n +2 > creh1.windows.txt
   fi
   python $GENOMICS_GENERAL/popgenWindows.py --windType predefined --windCoords creh1.windows.txt -f phased \
      -g erin63.geno.gz -o erin63.PopGenStats.csv --popsFile erin63.popmap.txt -p roumanicus -p europaeus -T 20 --writeFailedWindows
   rm creh1.windows.txt
fi

if [ ! -e windowStats.png ]; then
   R --slave --no-save < plot_windowStats.R
fi

# The previous estimation of population genetic parameters was done in the same windows where abbababa statistics were
# also estimated, for the sake of comparison. However, it is convenient to calculate population genetic statistics
# in non-overlapping windows, of a smaller size maybe, so that correlations among statistics can be evaluated
# without violating the assumption of independence among values.

if [ ! -e erin63.AccuStats.csv ]; then
   python $GENOMICS_GENERAL/popgenWindows.py -g erin63.geno.gz -o erin63.AccuStats.csv -f phased --windType sites \
      -w 50 -m 40 -O 0 --popsFile erin63.popmap.txt -p roumanicus -p europaeus -T 20 --writeFailedWindows
fi

# At the same time, I also want to be able to compare population genetics statistics and abba/baba statistics estimated
# in non-overlapping and common windows, even though from different subsets of SNPs.
if [ ! -e nonoverlap.abbababa.csv ]; then
   python $GENOMICS_GENERAL/ABBABABAwindows.py -g creh1.geno.gz -f phased -o nonoverlap.abbababa.csv --windType sites \
   -w 110 --overlap 0 -m 100 --popsFile creh1.popmap.txt -P1 concolor -P2 roumanicus -P3 europaeus -O Hemiechinus --writeFailedWindows -T 20
fi
if [ ! -e nonoverlap.PopGenStats.csv ]; then
   if [ ! -e nonoverlap.windows.txt ]; then
      cut -f 1,2,3 -d ',' --output-delimiter=$'\t' nonoverlap.abbababa.csv | tail -n +2 > nonoverlap.windows.txt
   fi
   python $GENOMICS_GENERAL/popgenWindows.py --windType predefined --windCoords nonoverlap.windows.txt -f phased \
   -g erin63.geno.gz -o nonoverlap.PopGenStats.csv --popsFile erin63.popmap.txt -p roumanicus -p europaeus -T 20 --writeFailedWindows
fi
if [ ! -e pi_dxy_fd.png ]; then
   R --slave --no-save < plot_pi_dxy_fd.R
   if [ -e nonoverlap.windows.txt ]; then rm nonoverlap.windows.txt; fi
fi
# Genes and intergenic regions are disjoint windows, but I can't separate them using regular windows, and I want
# to estimate parameters separately for the two kinds of regions. Thus, below I split the geno files. Warning:
# I assume that whatever does not overlap a gene is intergenic.

if [ ! -e ../../data/annotation.gff3.gz ]; then
   wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Erinaceus_europaeus/GFF/ref_EriEur2.0_top_level.gff3.gz
   mv ref_EriEur2.0_top_level.gff3.gz ../../data/annotation.gff3.gz
fi

if [ ! -e genes.bed ]; then
   gunzip -c ../../data/annotation.gff3.gz | grep -P "^#|\tgene\t" | bedtools merge -i - > genes.bed
fi

if [ ! -e genes.creh1.geno.gz ]; then
   if [ ! -e genic_filter.creh1.txt ]; then
      if [ ! -e snps.creh1.bed ]; then
         gunzip -c creh1.geno.gz | gawk '(NR == 1){print $0 >"genes.creh1.geno"}(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' > snps.creh1.bed
      fi
      bedtools intersect -a snps.creh1.bed -b genes.bed -wa | gawk '{print $1 "\t" $3 "\t"}' > genic_filter.creh1.txt
   fi
   gunzip -c creh1.geno.gz | grep -F -f genic_filter.creh1.txt >> genes.creh1.geno
   gzip genes.creh1.geno
fi

if [ ! -e inter.creh1.geno.gz ]; then
   if [ ! -e genic_filter.creh1.txt ]; then
      if [ ! -e snps.creh1.bed ]; then
         gunzip -c creh1.geno.gz | gawk '(NR == 1){print $0 >"inter.creh1.geno"}(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' > snps.creh1.bed
      fi
      bedtools intersect -a snps.creh1.bed -b genes.bed -wa | gawk '{print $1 "\t" $3 "\t"}' > genic_filter.creh1.txt
   fi
   gunzip -c creh1.geno.gz | grep -F -v -f genic_filter.creh1.txt >> inter.creh1.geno
   gzip inter.creh1.geno
fi

if [ ! -e genes.erin63.geno.gz ]; then
   gunzip -c erin63.geno.gz | head -n 1 > genes.erin63.geno
   if [ ! -e genic_filter.erin63.txt ]; then
      if [ ! -e snps.erin63.bed ]; then
         gunzip -c erin63.geno.gz | gawk '(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' > snps.erin63.bed
      fi
      bedtools intersect -a snps.erin63.bed -b genes.bed -wa | gawk '{print $1 "\t" $3 "\t"}' > genic_filter.erin63.txt
   fi
   gunzip -c erin63.geno.gz | grep -F -f genic_filter.erin63.txt >> genes.erin63.geno
   gzip genes.erin63.geno
fi

if [ ! -e inter.erin63.geno.gz ]; then
   gunzip -c erin63.geno.gz | head -n 1 > inter.erin63.geno
   if [ ! -e genic_filter.erin63.txt ]; then
      if [ ! -e snps.erin63.bed ]; then
         gunzip -c erin63.geno.gz | gawk '(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' > snps.erin63.bed
      fi
      bedtools intersect -a snps.erin63.bed -b genes.bed -wa | gawk '{print $1 "\t" $3 "\t"}' > genic_filter.erin63.txt
   fi
   gunzip -c erin63.geno.gz | grep -F -v -f genic_filter.erin63.txt >> inter.erin63.geno
   gzip inter.erin63.geno
   rm genic_filter.erin63.txt
fi

# Below I calculate abba/baba and population genetic statistics separately for genic and intergenic
# regions. Now, I am more interested in genome-wide values than in variation along contigs. It does
# not make much sense to use the same windows in both analyses, unless I used just the whole contigs.
# However, instead of running the estimates in differently-sized contigs, that then must be averaged,
# I prefer to use windows with a comparable amount of SNPs each, without overlap.

if [ ! -e genes.abbababa.csv ]; then
   python $GENOMICS_GENERAL/ABBABABAwindows.py -g genes.creh1.geno.gz -f phased -o genes.abbababa.csv --windType sites \
      -w 110 --overlap 0 -m 100 --popsFile creh1.popmap.txt -P1 concolor -P2 roumanicus -P3 europaeus -O Hemiechinus --writeFailedWindows -T 20
fi

if [ ! -e genes.PopGenStats.csv ]; then
   python $GENOMICS_GENERAL/popgenWindows.py --windType sites -w 110 --overlap 0 -m 100 -f phased \
      -g genes.erin63.geno.gz -o genes.PopGenStats.csv --popsFile erin63.popmap.txt -p roumanicus -p europaeus -T 20 --writeFailedWindows
fi

if [ ! -e inter.abbababa.csv ]; then
   python $GENOMICS_GENERAL/ABBABABAwindows.py -g inter.creh1.geno.gz -f phased -o inter.abbababa.csv --windType sites \
      -w 110 --overlap 0 -m 100 --popsFile creh1.popmap.txt -P1 concolor -P2 roumanicus -P3 europaeus -O Hemiechinus --writeFailedWindows -T 20
fi

if [ ! -e inter.PopGenStats.csv ]; then
   python $GENOMICS_GENERAL/popgenWindows.py --windType sites -w 110 --overlap 0 -m 100 -f phased \
      -g inter.erin63.geno.gz -o inter.PopGenStats.csv --popsFile erin63.popmap.txt -p roumanicus -p europaeus -T 20 --writeFailedWindows
fi

if [ ! -e pi_dxy.png ]; then
   R --slave --no-save < plot_pi_dxy.R
fi

if [ ! -e D.png ]; then
   R --slave --no-save < plot_D.R
fi

# I expected nucleotide diversity (pi) in either E. roumanicus or E. europaeus to be proportional to
# divergence between the two. The rationale is that under neutrality diversity is expected to be 4Nu
# and divergence 2Tu. Thus, divergence should be proportional to diversity, with a coefficient similar
# to T/2N, where T is divergence time and N, population size.
#
# What I see is that divergence is negatively correlated with diversity: regions with high diversity,
# more in E. roumanicus than in E. europaeus, are typically less diverged, more conserved. The effect
# is small, but significant. This does not make any sense from the neutral point of view. The pattern
# would be better explained by positive selection driving the removal of variation and the fixation
# of differences (divergence).
#
# Another possible explanation is introgression by itself, because it adds variation at the same time that
# it reduces the genetic distance between the species. The pattern also fits with the idea that introgression
# is polarized toward E. roumanicus.
#
# I do not plot fd, because it only makes sense for the windows where D > 0. Actually, instead of windows,
# it is clear that I should estimated an fd for all intergenic and one for all genic regions. The same strategy
# used in 2018-10-24 should be useful here. Note that to estimate f_d, I do not need to split the P3 population
# in two subpopulations. But, to re-use the 'estimate_fd.R' script, I should add dummy variables in place of P3a
# and P3b.

if [ ! -e D.txt ] || [ ! -e fd.txt ]; then
   if [ ! -e genic.tsv ] || [ ! -e inter.tsv ]; then
      if [ ! -e creh1.tsv ]; then
         gawk -v P1='concolor' \
              -v P2='roumanicus' \
              -v P3='europaeus' \
              -v OUTGROUP='Hemiechinus' \
              -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MINOUT=1 \
              -f $FREQ2_AWK creh1.popmap.txt $VCF1 | grep -v 999.9999 > creh1.tsv
      fi

      if [ ! -e genic_filter.txt ]; then
         if [ ! -e snps.bed ]; then
            gawk '(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' creh1.tsv > snps.bed
         fi
         bedtools intersect -a snps.bed -b genes.bed -wa | gawk '{print $1 "\t" $3 "\t"}' > genic_filter.txt
         #rm snps.bed
      fi

      if [ ! -e genic.tsv ]; then
         echo -e "scaffold\tposition\tconcolor\troumanicus\tdummy1\tdummy2\teuropaeus" > genic.tsv
         grep -F -f genic_filter.txt creh1.tsv | gawk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tNA\tNA\t" $5}' >> genic.tsv
      fi

      if [ ! -e inter.tsv ]; then
         echo -e "scaffold\tposition\tconcolor\troumanicus\tdummy1\tdummy2\teuropaeus" > inter.tsv
         grep -F -v -f genic_filter.txt creh1.tsv | gawk '(NR > 1){print $1 "\t" $2 "\t" $3 "\t" $4 "\tNA\tNA\t" $5}' >> inter.tsv
         #rm genic_filter
      fi
   fi
   #rm creh1.tsv

   if [ ! -e D.txt ]; then
      if [ ! -e contig_lengths.txt ]; then
         gawk '(/^##contig=<ID=/){split($1,A,/[=,>]/); print A[3] "\t" A[5]}' $VCF1 > contig_lengths.txt
      fi
      echo "# D statistic in genic regions:" > D.txt
      R -q --no-save <abba_baba.R --args genic.tsv contig_lengths.txt concolor roumanicus europaeus 1e5 | grep -vP "^[#>\+]" >> D.txt
      echo "# D statistic in intergenic regions:" >> D.txt
      R -q --no-save <abba_baba.R --args inter.tsv contig_lengths.txt concolor roumanicus europaeus 1e5 | grep -vP "^[#>\+]" >> D.txt
      rm contig_lengths.txt
   fi

   if [ ! -e fd.txt ]; then
      echo "# f_d statistic in genic regions:" > fd.txt
      R -q --no-save <estimate_f.R --args genic.tsv concolor roumanicus dummy1 dummy2 europaeus | grep -vP "^[#>\+]" >> fd.txt
      echo "# f_d statistic in intergenic regions:" >> fd.txt
      R -q --no-save <estimate_f.R --args inter.tsv concolor roumanicus dummy1 dummy2 europaeus | grep -vP "^[#>\+]" >> fd.txt
   fi
fi

#
# +------------+----------+-----------+-----------+-----------+
# |   Region   |     D    |   D_err   |   f_hom   |    fd     |
# +------------+----------+-----------+-----------+-----------+
# |      Genic | 0.108530 | 0.0278940 | 0.0264087 | 0.0157880 |
# | Intergenic | 0.170063 | 0.0230392 | 0.0436104 | 0.0255906 |
# +------------+----------+-----------+-----------+-----------+
#
# As expected, genic regions have less introgression than intergenic ones, confirming that
# introgression is in general deleterious.

if [ ! -e genes_in_extreme_fd_windows.gff ]; then
   if [ ! -e extreme_fd.bed ]; then
      # In nonoverlap.abbababa.csv, fields 10 and 11 are fd and fdM. I determined manually the thresholds
      # 0.20 and 0.15, respectively, to select the windows with the largest amount of introgression.
      gawk -v FS=',' '((NR > 1) && ($10 !~ /nan/) && (($10 > 0.20) || ($11 > 0.15))){
         print $1 "\t" $2 - 1 "\t" $3
      }' nonoverlap.abbababa.csv > extreme_fd.bed
   fi
   gunzip -c ../../data/annotation.gff3.gz | \
   grep -P "\tgene\t" | \
   bedtools intersect -a extreme_fd.bed -b - -F 1.0 -wb | grep -v pseudogen > genes_in_extreme_fd_windows.gff
fi

# There are about 100 protein coding genes in the windows with the highest proportions of
# introgressed variation. Because natural selection seems to be removing introgressed variation,
# we expect that if those genes harbor any introgressed variatio (it could be around them instead),
# then they must be the most tolerant ones to introgression. The point is, we don't have reasons
# to believe that any of those genes contain 'adaptive' introgression. Thus, a functional enrichment
# analysis may be unnecessary: I do not expect any functional category to be more tolerant to
# maladaptive introgression than any other. And if there was one, it would not be very interesting,
# anyways.
#
# What would be interesting is to find signatures of positive selection in any gene located in a
# highly introgressed region. After having informally rejected neutrality in E. roumanicus through
# a genome-wide version of the Hudson-Kreitman-Aguadé test, now we could apply a McDonald-Kreitman
# test to every gene. For that, we need tables of synonymous and non-synonymous polymorphisms and
# fixed differences betwen the two species. We can obtain that information from the vcf file and from
# the populations map file.
#
# The McDonald-Kreitman test would also be very informative about one window with the highest level
# of divergence, dxy=0.6172, which happens to have one of the lowest levels of diversity within E.
# europaeus (0.0523, in nonoverlap.PopGenStats.csv). The file erin63.AccuStats.csv has two windows
# in the same scaffold, NW_006804247.1, with similar values. The whole scaffold could have experienced
# a recent selective sweep.
#
# I tried the program SnpEff, to annotate the functional effects of SNPs in the vcf. However, it uses
# pre-configured annotation files from Ensembl, which are based on a different assembly. The scaffold
# names are different, and it is not compatible with our dataset.
