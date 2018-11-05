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

# Different vcf files are optimized for different statistics:
VCF1=../2018-09-26/ErinaceusAndHemiechinus_r1e1c1h1.vcf
VCF2=../2018-09-26/ErinMaxMiss63_r10e10c4.vcf
# This is the folder where Simon H. Martin's genomic_general repository is cloned:
GENOMICS_GENERAL=~/bin/genomics_general
POPDATA=../../data/populations.txt

# I will exclude both the hybrid and the other admixed individual, in order to work only
# with the historical signal of introgression. It will help interpretation.

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
      gunzip -c erin63.geno.gz | head -n 1 | sed 's/\t/\n/g' | gail -n +3 > good_samples.txt
      grep -F -f good_samples.txt popmap.txt > erin63.popmap.txt
   fi
fi

if [ ! -e erin63.PopGenStats.csv ]; then
   python $GENOMICS_GENERAL/popgenWindows.py --windType sites -w 50 -m 50 -O 40 -f phased \
      -g erin63.geno.gz -o erin63.PopGenStats.csv --popsFile erin63.popmap.txt -p roumanicus -p europaeus -T 20
fi

if [ ! -e creh1.abbababa.csv ]; then

fi
