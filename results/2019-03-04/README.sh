#!/bin/bash
#
#				2019-03-04
#				----------
#
# On 2019-03-01, I obtained tables of polymorphism and divergence to run the asymptotic
# McDonald-Kreitman test. However, I did not identify the derived allele. That was necessary
# to distribute the polymorphic sites amont derived-allele frequency classes. To identify
# the derived allele, I need to include the information of the Hemiechinus sample, and
# assume that the allele present there is the ancestral one.

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
   echo Er73_SNG1      > undesired.txt    # Atelerix
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
#
# What matters is to estimate allele frequencies in the focal species, where I will request
# a minimum number of individuals. In the sister species, one individual should be enough to
# determine if the site is divergent (i.e., all individuals sampled in the focal species homozygous
# for an allele different to that present in the sister species' indiviual(s)).
for i in 7 8 9 10 11 12; do
   if [ ! -e R${i}E1.polymorphism.txt ] || [ ! -e R${i}E1.divergence.txt ] || [ ! -e E${i}R1.polymorphism.txt ] || [ ! -e E${i}R1.divergence.txt ]; then
      if [ ! -e R${i}E1.vcf ] || [ ! -e E${i}R1.vcf ]; then
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
         if [ ! -e R${i}E1.vcf ]; then
            gawk -v ROU=$i -v EUR=1 -v CON=0 -v HEM=1 -f ../../bin/filtervcf.awk popmap.txt flagged.vcf > R${i}E1.vcf
         fi
         if [ ! -e E${i}R1.vcf ]; then
            gawk -v ROU=1 -v EUR=$i -v CON=0 -v HEM=1 -f ../../bin/filtervcf.awk popmap.txt flagged.vcf > E${i}R1.vcf
         fi
      fi
      if [ ! -e R${i}E1.polymorphism.txt ] || [ ! -e R${i}E1.divergence.txt ]; then
         gawk -f ../../bin/vcf2MK_2.awk -v FOCAL="roumanicus" -v SISTER="europaeus" -v OUTPUT1="R${i}E1.polymorphism.txt" -v OUTPUT2="R${i}E1.divergence.txt" popmap.txt R${i}E1.vcf
      fi
      if [ ! -e E${i}R1.polymorphism.txt ] || [ ! -e E${i}R1.divergence.txt ]; then
         gawk -f ../../bin/vcf2MK_2.awk -v FOCAL="europaeus" -v SISTER="roumanicus" -v OUTPUT1="E${i}R1.polymorphism.txt" -v OUTPUT2="E${i}R1.divergence.txt" popmap.txt E${i}R1.vcf
      fi
   fi
done
if [ -e flagged.vcf ]; then rm flagged.vcf; fi

if [ ! -e summary.txt ]; then
   echo -e "E. roumanicus" >  summary.txt
   echo -e "-------------" >> summary.txt
   echo -e " \tPolymorphism\tPolymorphism\tDivergence\tDivergence" >> summary.txt
   echo -e "N\tSynonymous\tNon-synonymous\tSynonymous\tNon-synonuymous\tAlpha" >> summary.txt
   for i in 7 8 9 10 11 12; do
      NON_POL=$(gawk '{S += $2}END{print S}' R${i}E1.polymorphism.txt)
      SYN_POL=$(gawk '{S += $3}END{print S}' R${i}E1.polymorphism.txt)
      SYN_DIV=$(head -n 1 R${i}E1.divergence.txt | cut -d " " -f 4)
      NON_DIV=$(tail -n 1 R${i}E1.divergence.txt | cut -d " " -f 4)
      ALPHA=$(echo "1 - ($SYN_DIV / $NON_DIV) * ($NON_POL / $SYN_POL)" | bc -l)
      echo -e "$i\t$SYN_POL      \t$NON_POL      \t$SYN_DIV     \t$NON_DIV     \t$ALPHA" >> summary.txt
   done
   echo >> summary.txt
   echo "E. europaeus" >> summary.txt
   echo "------------" >> summary.txt
   echo -e " \tPolymorphism\tPolymorphism\tDivergence\tDivergence" >> summary.txt
   echo -e "N\tSynonymous\tNon-synonymous\tSynonymous\tNon-synonuymous\tAlpha" >> summary.txt
   for i in 7 8 9 10 11 12; do
      NON_POL=$(gawk '{S += $2}END{print S}' E${i}R1.polymorphism.txt)
      SYN_POL=$(gawk '{S += $3}END{print S}' E${i}R1.polymorphism.txt)
      SYN_DIV=$(head -n 1 E${i}R1.divergence.txt | cut -d " " -f 4)
      NON_DIV=$(tail -n 1 E${i}R1.divergence.txt | cut -d " " -f 4)
      ALPHA=$(echo "1 - ($SYN_DIV / $NON_DIV) * ($NON_POL / $SYN_POL)" | bc -l)
      echo -e "$i\t$SYN_POL      \t$NON_POL      \t$SYN_DIV     \t$NON_DIV     \t$ALPHA" >> summary.txt
   done

fi

if [ ! -e asymptoticMK.html ]; then
   R -q --save -e "rmarkdown::render('asymptoticMK.Rmd', output_file='asymptoticMK.html')"
fi
