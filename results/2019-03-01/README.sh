#!/bin/bash
#
#				2019-03-01
#				----------
#
# On 2019-02-27 I annotated the vcf file. Now, the goal is to filter the annotated
# vcf and prepare the input file for an MK test. I aim at using the asymptotic MK
# test in http://benhaller.com/messerlab/asymptoticMK.html. I want to test neutrality
# in E. roumanicus and E. europaeus. The divergence will be between these two species.
# I will use only SNPs annotated as either 'synonymous' or 'missense'. I will filter
# the vcf file to include sites with at least 10 E. europaeus and 10 E. roumanicus
# samples.

VCF=../2019-02-27/merged.annotated.vcf
POPDATA=../../data/populations.txt

if [ ! -e popmap.txt ]; then
   gawk 'BEGIN{
      POPNAME[1] = "roumanicus"
      POPNAME[2] = "europaeus"
      POPNAME[3] = "concolor"
      POPNAME[4] = "hybrid"
      POPNAME[5] = "Hemiechinus"
      POPNAME[6] = "Atelerix"
   }(/^[^#]/){
      print $1 "\t" POPNAME[$2]
   }' $POPDATA | sort -k 2,2 > popmap.txt
fi

if [ ! -e undesired.txt ]; then
   echo Er65_IS25      > undesired.txt    # Hemiechinus
   echo Er73_SNG1     >> undesired.txt    # Atelerix
   echo Er59_FR1      >> undesired.txt    # wrong
   echo Er68_PRT1B    >> undesired.txt    # wrong
   echo Er38_LB1      >> undesired.txt    # concolor
   echo Er49_GR38     >> undesired.txt    # concolor
   echo Er53_ASR7     >> undesired.txt    # concolor
   echo Er64_IS1      >> undesired.txt    # concolor
   echo Er75_TRC2A    >> undesired.txt    # concolor
   echo Er27_SK32     >> undesired.txt    # admixed concolor
   echo Er26_JUG4     >> undesired.txt    # admixed concolor
fi

# This is a sites-kind of analysis, rather than an individual-kind of analysis. That is,
# I don't need individuals to have complete data, and there is no need to exclude individuals.
# Thus, I could add the binary presence flag to the original vcf, and then filter sites.
# However, the input file is large and dirty. It may even be faster (if not only cleaner)
# to remove sites with indels or with more than one allele, and to exclude individuals
# not targeted in this analysis before adding the flag.

for i in 7 8 9 10 11 12; do
   if [ ! -e R${i}E${i}.polymorphism.txt ] || [ ! -e R${i}E${i}.divergence.txt ]; then
      if [ ! -e R${i}E${i}.vcf ]; then
         if [ ! -e flagged.vcf ] ; then
            vcftools --vcf $VCF \
                     --stdout \
                     --remove undesired.txt \
                     --maf 0.01 \
                     --max-maf 0.95 \
                     --remove-indels \
                     --min-alleles 2 \
                     --max-alleles 2 \
                     --maxDP 200 \
                     --minDP 4 \
                     --minQ 50 \
                     --thin 261 \
                     --recode \
                     --recode-INFO-all | \
            gawk -f ../../bin/add_flag.awk > flagged.vcf
         fi
         gawk -v ROU=10 -v EUR=10 -v CON=0 -f ../../bin/filtervcf.awk popmap.txt flagged.vcf > R${i}E${i}.vcf
      fi
      gawk -f ../../bin/vcf2MK.awk -v FOCAL="roumanicus" -v SISTER="europaeus" -v OUTPUT1="R${i}E${i}.polymorphism.txt" -v OUTPUT2="R${i}E${i}.divergence.txt" popmap.txt R${i}E${i}.vcf
      gawk -f ../../bin/vcf2MK.awk -v FOCAL="europaeus" -v SISTER="roumanicus" -v OUTPUT1="E${i}R${i}.polymorphism.txt" -v OUTPUT2="E${i}R${i}.divergence.txt" popmap.txt R${i}E${i}.vcf
      # rm R${i}E${i}.vcf
   fi
done
rm flagged.vcf
