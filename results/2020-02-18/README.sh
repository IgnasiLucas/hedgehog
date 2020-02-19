#!/bin/bash
#
#				2020-02-18
#				----------
#
# On 2019-03-04 I estimated the proportion of fixed differencese between E. roumanicus
# and E. europaeus driven by positive selection, using script vcf2mk_2.awk to count the
# number of sites of each category. Kirstýna reviewed those results and they are different
# from others. I want to make sure that the scripts worked as expected.

VCF=../2019-02-27/merged.annotated.vcf
POPDATA=../../data/populations.txt
MINSAMPLE=7   # Set here the minimum number of individuals per population

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

# Below I filter the original vcf file before adding the BPF. Then, I filter again to
# retain only sites where at least MINSAMPLE E. roumanicus, MINSAMPLE E. europaeus, and 1 Hemiechinus
# individuals are genotyped. Note that the exact number of individuals genotyped in a
# locus varies, and that's why I have to measure derived allele frequency with a float,
# instead of using the integer number of derived allele counts (see vcf2MK_2.awk).

if [ ! -d original ]; then mkdir original; fi
for i in $MINSAMPLE; do
   if [ ! -e original/R${i}E${i}.polymorphism.txt ] || [ ! -e original/R${i}E${i}.divergence.txt ]; then
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
         gawk -v ROU=$i -v EUR=$i -v CON=0 -v HEM=1 -f ../../bin/filtervcf.awk popmap.txt flagged.vcf > R${i}E${i}.vcf
      fi
      gawk -f vcf2MK_3.awk -v FOCAL="roumanicus" -v SISTER="europaeus" -v OUTPUT1="original/R${i}E${i}.polymorphism.txt" -v OUTPUT2="original/R${i}E${i}.divergence.txt" -v OUTPUT3="original/frequencies_rou.txt" popmap.txt R${i}E${i}.vcf
      gawk -f vcf2MK_3.awk -v FOCAL="europaeus" -v SISTER="roumanicus" -v OUTPUT1="original/E${i}R${i}.polymorphism.txt" -v OUTPUT2="original/E${i}R${i}.divergence.txt" -v OUTPUT3="original/frequencies_eur.txt" popmap.txt R${i}E${i}.vcf
   fi
done
if [ -e flagged.vcf ]; then rm flagged.vcf; fi

if [ ! -e original/summary.txt ]; then
   echo -e "E. roumanicus" >  original/summary.txt
   echo -e "-------------" >> original/summary.txt
   echo -e " \tPolymorphism\tPolymorphism\tDivergence\tDivergence" >> original/summary.txt
   echo -e "N\tSynonymous\tNon-synonymous\tSynonymous\tNon-synonuymous\tAlpha" >> original/summary.txt
   for i in $MINSAMPLE; do
      NON_POL=$(gawk '{S += $2}END{print S}' original/R${i}E${i}.polymorphism.txt)
      SYN_POL=$(gawk '{S += $3}END{print S}' original/R${i}E${i}.polymorphism.txt)
      SYN_DIV=$(head -n 1 original/R${i}E${i}.divergence.txt | cut -d " " -f 4)
      NON_DIV=$(tail -n 1 original/R${i}E${i}.divergence.txt | cut -d " " -f 4)
      ALPHA=$(echo "1 - ($SYN_DIV / $NON_DIV) * ($NON_POL / $SYN_POL)" | bc -l)
      echo -e "$i\t$SYN_POL      \t$NON_POL      \t$SYN_DIV     \t$NON_DIV     \t$ALPHA" >> original/summary.txt
   done
   echo >> original/summary.txt
   echo "E. europaeus" >> original/summary.txt
   echo "------------" >> original/summary.txt
   echo -e " \tPolymorphism\tPolymorphism\tDivergence\tDivergence" >> original/summary.txt
   echo -e "N\tSynonymous\tNon-synonymous\tSynonymous\tNon-synonuymous\tAlpha" >> original/summary.txt
   for i in $MINSAMPLE; do
      NON_POL=$(gawk '{S += $2}END{print S}' original/E${i}R${i}.polymorphism.txt)
      SYN_POL=$(gawk '{S += $3}END{print S}' original/E${i}R${i}.polymorphism.txt)
      SYN_DIV=$(head -n 1 original/E${i}R${i}.divergence.txt | cut -d " " -f 4)
      NON_DIV=$(tail -n 1 original/E${i}R${i}.divergence.txt | cut -d " " -f 4)
      ALPHA=$(echo "1 - ($SYN_DIV / $NON_DIV) * ($NON_POL / $SYN_POL)" | bc -l)
      echo -e "$i\t$SYN_POL      \t$NON_POL      \t$SYN_DIV     \t$NON_DIV     \t$ALPHA" >> original/summary.txt
   done

fi

# Up to here, it's just what I did before. Now, I want to count sites in a different way, to see
# if I get the same result. I will start from the R10E10.vcf file, so that I will only be checking
# that vcf2MK_2.awk works as expected, which I think is the most obscure piece of code.

if [ ! -d new ]; then mkdir new; fi
if [ ! -e new/counts_table.txt ]; then
   gawk '(FILENAME ~ /popmap/){
      POPULATION[$1] = $2
   }((FILENAME ~ /vcf/) && (/^#CHROM/)){
      for (i=0; i<=(NF-10); i++){
         POP[i] = POPULATION[$(i+10)]
      }
      print("CHROM\tPOS\tEFFECT\troumanicus_0\troumanicus_1\troumanicus_2\teuropaeus_0\teuropaeus_1\teuropaeus_2\tHemiechinus_0\tHemiechinus_1\tHemiechinus_2")
   }((FILENAME ~ /vcf/) && (/synonymous_variant/) && (!(/missense_variant/))){
      delete(COUNT)
      for (i=0; i<=(NF-10); i++){
         split($(i+10), GT, /:/)
         COUNT[GT[1], POP[i]]++
      }
      print($1 "\t" $2 "\tsynonymous\t" COUNT["0/0", "roumanicus"]  + 0 "\t" COUNT["0/1", "roumanicus"]  + 0 "\t" COUNT["1/1", "roumanicus"]  + 0 "\t" \
                                        COUNT["0/0", "europaeus"]   + 0 "\t" COUNT["0/1", "europaeus"]   + 0 "\t" COUNT["1/1", "europaeus"]   + 0 "\t" \
                                        COUNT["0/0", "Hemiechinus"] + 0 "\t" COUNT["0/1", "Hemiechinus"] + 0 "\t" COUNT["1/1", "Hemiechinus"] + 0)
   }((FILENAME ~ /vcf/) && (/missense_variant/)){
      delete(COUNT)
      for (i=0; i<=(NF-10); i++){
         split($(i+10), GT, /:/)
         COUNT[GT[1], POP[i]]++
      }
      print($1 "\t" $2 "\tmissense\t" COUNT["0/0", "roumanicus"]  + 0 "\t" COUNT["0/1", "roumanicus"]  + 0 "\t" COUNT["1/1", "roumanicus"]  + 0 "\t" \
                                      COUNT["0/0", "europaeus"]   + 0 "\t" COUNT["0/1", "europaeus"]   + 0 "\t" COUNT["1/1", "europaeus"]   + 0 "\t" \
                                      COUNT["0/0", "Hemiechinus"] + 0 "\t" COUNT["0/1", "Hemiechinus"] + 0 "\t" COUNT["1/1", "Hemiechinus"] + 0)
   }' popmap.txt R${MINSAMPLE}E${MINSAMPLE}.vcf > new/counts_table.txt
fi

# Table 'counts_table.txt' contains all the information needed to estimate the diverngece
# and polymorphism levels in both missense and synonymous substitutions. This intermediate
# file makes it easier to check the numbers manually. We need the following information:
#
#  - d_0: divergence at synonymous sites. I count divergent sites as those where different
#         alleles are fixed in the focal and sister populations. Sites where the focal
#         population is fixed and the sister one is polymorphic contribute partially to
#         divergence, to the extent of the frequency of the non-focal allele in the sister.
#  - d:   divergence at missense sites. Similarly.
#  - p_0: number of polymorphic synonymous sites with frequency x.
#  - p:   number of polymorphic missense sites with frequency x.
#  - x:   If genotypes where available for all samples in all sites, x could take all
#         possible values of the number of alleles seen in the population. But I am accepting
#         sites with a different number of individuals genotyped. Thus, I estimate allele
#         frequency as a real number (a float) in every site, and I bin x by rounding. This
#         could introduce a bias, and it would be reasonable to use the bins that could
#         be observed in the minimum number of individuals accepted. For example, if we
#         use sites with down to 7 individuals (14 alleles), possible values with the
#         minimum number of individuals are 1/14, 2/14, 3/14... 14/14.

if [ ! -e new/frequencies.txt ]; then
   # I exclude sites where the outgroup is heterozygous: unclear what's the ancestral allele
   gawk '((NR > 1) && ($11 + 0 == 0)){
      if ($10 + 0 == 1) {
         # The derived allele is the alternative
         ROUMANICUS = ($5 + 2.0 * $6) / (2.0 * ($4 + $5 + $6))   # derived allele frequency in E. roumanicus
         EUROPAEUS  = ($8 + 2.0 * $9) / (2.0 * ($7 + $8 + $9))   # derived allele frequency in E. europaeus
      } else {
         # The derived allele is the reference
         ROUMANICUS = ($5 + 2.0 * $4) / (2.0 * ($4 + $5 + $6))
         EUROPAEUS  = ($8 + 2.0 * $7) / (2.0 * ($7 + $8 + $9))
      }
      print $1 "\t" $2 "\t" $3 "\t" ROUMANICUS "\t" EUROPAEUS
   }' new/counts_table.txt > new/frequencies.txt
fi

# The frequencies.txt table shows the real frequency estimates, and allows
# us to customize the binning of frequencies. For example, I can bin frequencies
# in 8 classes: (0,0.125], (0.125,0.250], (0.250,0.375], (0.375,0.500],
# (0.500,0.625], (0.625,0.750], (0.750,0.875], (0.875,1.000):

if [ ! -e new/R${MINSAMPLE}E${MINSAMPLE}.polymorphism.txt ]; then
   gawk -v OUTPUT1="new/R${MINSAMPLE}E${MINSAMPLE}.polymorphism.txt" -v OUTPUT2="new/R${MINSAMPLE}E${MINSAMPLE}.divergence.txt" '{
      if ($4 + 0 == 0.0) {
         DIVERGENCE[$3] += $5
      } else {
         if ($4 + 0 < 0.125) {
            POLYMORPHISM[$3, "0.125"]++
         } else {
            if ($4 + 0 < 0.250) {
               POLYMORPHISM[$3, "0.250"]++
            } else {
               if ($4 + 0 < 0.375) {
                  POLYMORPHISM[$3, "0.375"]++
               } else {
                  if ($4 + 0 < 0.500) {
                     POLYMORPHISM[$3, "0.500"]++
                  } else {
                     if ($4 + 0 < 0.625) {
                        POLYMORPHISM[$3, "0.625"]++
                     } else {
                        if ($4 + 0 < 0.750) {
                           POLYMORPHISM[$3, "0.750"]++
                        } else {
                           if ($4 + 0 < 0.875) {
                              POLYMORPHISM[$3, "0.875"]++
                           } else {
                              if ($4 + 0 < 1.000) {
                                 POLYMORPHISM[$3, "1.000"]++
                              } else {
                                 DIVERGENCE[$3] += 1 - $5
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   } END{
#     print "Frequency\tMissense\tSynonymous"
      for (i = 0.125; i <= 1.000; i += 0.125) {
#        printf("%.3f-%.3f\t%u\t%u\n", i-0.125, i, POLYMORPHISM["missense", sprintf("%.3f",i)], POLYMORPHISM["synonymous", sprintf("%.3f", i)])
         printf("%.4f\t%u\t%u\n", i - 0.0625, POLYMORPHISM["missense", sprintf("%.3f",i)], POLYMORPHISM["synonymous", sprintf("%.3f", i)]) >> OUTPUT1
      }
      print "Missense divergence: "   DIVERGENCE["missense"]    > OUTPUT2
      print "Synonymous divergence: " DIVERGENCE["synonymous"] >> OUTPUT2
   }' new/frequencies.txt
fi


if [ ! -e new/E${MINSAMPLE}R${MINSAMPLE}.polymorphism.txt ]; then
   gawk -v OUTPUT1="new/E${MINSAMPLE}R${MINSAMPLE}.polymorphism.txt" -v OUTPUT2="new/E${MINSAMPLE}R${MINSAMPLE}.divergence.txt" '{
      if ($5 + 0 == 0.0) {
         DIVERGENCE[$3] += $4
      } else {
         if ($5 + 0 < 0.125) {
            POLYMORPHISM[$3, "0.125"]++
         } else {
            if ($5 + 0 < 0.250) {
               POLYMORPHISM[$3, "0.250"]++
            } else {
               if ($5 + 0 < 0.375) {
                  POLYMORPHISM[$3, "0.375"]++
               } else {
                  if ($5 + 0 < 0.500) {
                     POLYMORPHISM[$3, "0.500"]++
                  } else {
                     if ($5 + 0 < 0.625) {
                        POLYMORPHISM[$3, "0.625"]++
                     } else {
                        if ($5 + 0 < 0.750) {
                           POLYMORPHISM[$3, "0.750"]++
                        } else {
                           if ($5 + 0 < 0.875) {
                              POLYMORPHISM[$3, "0.875"]++
                           } else {
                              if ($5 + 0 < 1.000) {
                                 POLYMORPHISM[$3, "1.000"]++
                              } else {
                                 DIVERGENCE[$3] += 1 - $4
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   } END{
#     print "Frequency\tMissense\tSynonymous"
      for (i = 0.125; i <= 1.000; i += 0.125) {
#        printf("%.3f-%.3f\t%u\t%u\n", i-0.125, i, POLYMORPHISM["missense", sprintf("%.3f",i)], POLYMORPHISM["synonymous", sprintf("%.3f", i)])
         printf("%.4f\t%u\t%u\n", i - 0.0625, POLYMORPHISM["missense", sprintf("%.3f",i)], POLYMORPHISM["synonymous", sprintf("%.3f", i)]) >> OUTPUT1
      }
      print "Missense divergence: "   DIVERGENCE["missense"]   >  OUTPUT2
      print "Synonymous divergence: " DIVERGENCE["synonymous"] >> OUTPUT2
   }' new/frequencies.txt
fi

# CONCLUSION
# ----------
#
# There was a small mistake in the vcf2MK_2.awk script. It was not reporting any synonymous polymorphisms at the
# lowest nor at the highest frequency classes (rounded as 0 and 1, respectively). After correcting for that, the
# values of alpha come closer to 0, but still negative.
#
# In the 'new' folder, I use a different binning of frequencies for the assymptotic estimation of alpha.
