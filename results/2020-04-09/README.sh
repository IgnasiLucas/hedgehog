#!/bin/bash
#
#				2020-04-09
#				==========
#
# The reference-genome bias works through a reduction in coverage in more
# divergent regions, which would downwardly bias both divergence and diversity
# measures, potentially introducing a negative correlation between them. My
# strategy here is to use sample-specific coverage data to assess and subtract
# the effect of coverage from the parameter estimates, using a linear model.
#
# Reproducibility notes: use an environment with python 2 and numpy.

# -----------------------------------------------------------------------------
#				  SET UP
# -----------------------------------------------------------------------------

VCF=../2018-09-26/ErinMaxMiss63_r10e10c4.vcf
GG=~/bin/genomics_general
POPDATA=../../data/populations.txt

# -----------------------------------------------------------------------------
#				   DATA
# -----------------------------------------------------------------------------

if [ ! -e erin63.popmap.txt ]; then
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

   if [ ! -e erin63.popmap.txt ]; then
      grep "^#CHROM" $VCF | sed 's/\t/\n/g' | tail -n +10 > good_samples.txt
      grep -F -f good_samples.txt popmap.txt > erin63.popmap.txt
      rm good_samples.txt
      rm popmap.txt
   fi
fi

for pop in concolor roumanicus europaeus; do
   if [ ! -e $pop.txt ]; then
      grep $pop erin63.popmap.txt | cut -f 1 > $pop.txt
   fi
   if [ ! -e $pop.ldepth ]; then
      vcftools --vcf $VCF --out $pop --keep $pop.txt --site-depth
   fi
   if [ ! -e $pop.frq ]; then
      vcftools --vcf $VCF --out $pop --keep $pop.txt --freq2
   fi
   if [ ! -e $pop.sites.pi ]; then
      vcftools --vcf $VCF --out $pop --keep $pop.txt --site-pi
   fi
done

if [ ! -e concolor_roumanicus.weir.fst ]; then
   vcftools --vcf $VCF --out concolor_roumanicus --keep concolor.txt --keep roumanicus.txt --weir-fst-pop concolor.txt --weir-fst-pop roumanicus.txt
fi

if [ ! -e concolor_europaeus.weir.fst ]; then
   vcftools --vcf $VCF --out concolor_europaeus --keep concolor.txt --keep europaeus.txt --weir-fst-pop concolor.txt --weir-fst-pop europaeus.txt
fi

if [ ! -e roumanicus_europaeus.weir.fst ]; then
   vcftools --vcf $VCF --out roumanicus_europaeus --keep roumanicus.txt --keep europaeus.txt --weir-fst-pop roumanicus.txt --weir-fst-pop europaeus.txt
fi

if [ ! -e data.tsv ]; then
   echo -e "SCAF\tPOS\tN_con\tN_rou\tN_eur\tAF_con\tAF_rou\tAF_eur\tDP_con\tDP_rou\tDP_eur\tPI_con\tPI_rou\tPI_eur\tFST_con_rou\tFST_con_eur\tFST_rou_eur" > data.tsv
   paste concolor.frq roumanicus.frq europaeus.frq 
fi

# -----------------------------------------------------------------------------
# 			REPORT
# -----------------------------------------------------------------------------

if [ ! -e reference.pdf ]; then
   R -q --no-save -e "rmarkdown::render('reference.Rmd', output_file='reference.pdf')
fi
