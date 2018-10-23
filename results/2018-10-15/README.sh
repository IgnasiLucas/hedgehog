#!/bin/bash
#
#				2018-10-15
#				==========
#
# Here I want to run the abba/baba or D test to compare the introgression signal
# among E. europaeus and E. roumanicus with the regions of the hybrid individual
# known to be admixed.

ADMIXTURE_DIR=/data/kristyna/hedgehog/results_2018/05-06-2018/r10e10c4
VCF=../2018-09-26/ErinaceusAndHemiechinus_r1e1c1h1.vcf
POPDATA=../../data/populations.txt
FREQ2_AWK=../../bin/freq2.awk

if [ ! -e popmap.txt ]; then
   # This file will still include the problematic individuals, but it
   # should not matter.
   gawk 'BEGIN{
      POPNAME[1] = "roumanicus"
      POPNAME[2] = "europaeus"
      POPNAME[3] = "concolor"
      POPNAME[4] = "hybrid"
      POPNAME[5] = "Hemiechinus"
      POPNAME[6] = "Atelerix"
   }(/^[^#]/){
      print $1 "\t" POPNAME[$2]
   }' $POPDATA > popmap.txt
fi

if [ ! -e popmap_split.txt ]; then
   # I will split E. europaeus according to the population structure
   # detected with Admixture. That means, there are only 4 E. europaeus
   # from the western subpopulation: Er30_79, Er74_SP16, Er63_IR6, and
   # Er56_AZ5. The four columns in Q4.txt represent ancestries from:
   # Eastern E. europaeus, E. concolor, E. roumanicus, and Western E.
   # europaeus, respectively. Note, that not all individuals in the current
   # vcf file were used for the Admixture analysis, and I need to add
   # them at the end.
   if [ ! -e Q4.txt ]; then
      paste $ADMIXTURE_DIR/out.012.indv $ADMIXTURE_DIR/erinaceus_41_r10e10c4.4.Q > Q4.txt
   fi
   gawk 'BEGIN{
      POP[2] = "western"
      POP[3] = "concolor"
      POP[4] = "roumanicus"
      POP[5] = "eastern"
   }{
      POPULATION = "none"
      for (i = 2; i <= 5; i++) {
         if ($i > 0.8) POPULATION = POP[i]
      }
      if (POPULATION == "none") POPULATION = "hybrid"
      print $1 "\t" POPULATION
   }END{
      print "Er36_SK24\troumanicus"
      print "Er58_FI7\teastern"      # This is a guess.
      print "Er61_GR87\troumanicus"
      print "Er62_GR95\troumanicus"
      print "Er67_IT5\teastern"      # This is a guess.
      print "Er65_IS25\tHemiechinus"
   }' Q4.txt > popmap_split.txt
fi

if [ ! -e creh1.tsv ]; then
   # "creh1" stands for "concolor, roumanicus, europaeus, Hemiechinus, 1 of each at least".
   gawk -v P1="concolor" \
        -v P2="roumanicus" \
        -v P3="europaeus" \
        -v OUTGROUP="Hemiechinus" \
        -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MINOUT=1 \
        -f $FREQ2_AWK popmap.txt $VCF > creh1.tsv
fi

if [ ! -e crweh1.tsv ]; then
   # "crweh1" stands for "concolor, roumanicus, western, eastern, Hemiechinus, 1 of each".
   gawk -v P1="concolor" \
        -v P2="roumanicus" \
        -v P3="eastern" \
        -v P4="western" \
        -v OUTGROUP="Hemiechinus" \
        -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MIN4=1 -v MINOUT=1 \
        -f $FREQ2_AWK popmap_split.txt $VCF | grep -v 999.9999 > crweh1.tsv
fi

# I found that Simon Martin estimates several flavours of the f statistic. For some of them
# we need the third population both split and not split. Thus, I shall combine creh1.tsv and
# crweh1.tsv in one file with the sites where all have data. The file crweh1.csv has fewer
# sites because both western and eastern populations are required to have data. All sites in
# crweh1.tsv are also present in creh1.tsv.
if [ ! -e combined.tsv ]; then
   cut -f 1,2 crweh1.tsv > zPositions.txt
   grep -F -f zPositions.txt creh1.tsv > z_creh1_subset.tsv
   paste crweh1.tsv z_creh1_subset.tsv | gawk '(($1 == $8) && ($2 == $9)){print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $12}' > combined.tsv
   rm zPositions.txt z_creh1_subset.tsv
fi

if [ ! -e contig_lengths.txt ]; then
   gawk '(/^##contig=<ID=/){split($1,A,/[=,>]/); print A[3] "\t" A[5]}' $VCF > contig_lengths.txt
fi

if [ ! -e D.txt ]; then
   echo -e "#Population_1\tPopulation_2\tPopulation_3\tD_statistic\tStandard_Error\tP-value\tNumber_Sites" > D.txt
fi
if ! grep -q -P "concolor\troumanicus\teuropaeus\t" D.txt; then
   R -q --no-save <abba_baba.R --args creh1.tsv contig_lengths.txt concolor roumanicus europaeus 1e5 | grep -vP "^[#>\+]" >> D.txt
fi
if ! grep -q -P "concolor\troumanicus\teastern\t" D.txt; then
   R -q --no-save <abba_baba.R --args crweh1.tsv contig_lengths.txt concolor roumanicus eastern 1e5 | grep -vP "^[#>\+]" >> D.txt
fi
if ! grep -q -P "concolor\troumanicus\twestern\t" D.txt; then
   R -q --no-save <abba_baba.R --args crweh1.tsv contig_lengths.txt concolor roumanicus western 1e5 | grep -vP "^[#>\+]" >> D.txt
fi
if ! grep -q -P "concolor\teastern\twestern\t" D.txt; then
   R -q --no-save <abba_baba.R --args crweh1.tsv contig_lengths.txt concolor eastern western 1e5 | grep -vP "^[#>\+]" >> D.txt
fi
if ! grep -q -P "western\teastern\tconcolor\t" D.txt; then
   R -q --no-save <abba_baba.R --args crweh1.tsv contig_lengths.txt western eastern concolor 1e5 | grep -vP "^[#>\+]" >> D.txt
fi
if ! grep -q -P "western\teastern\troumanicus\t" D.txt; then
   R -q --no-save <abba_baba.R --args crweh1.tsv contig_lengths.txt western eastern roumanicus 1e5 | grep -vP "^[#>\+]" >> D.txt
fi

if [ ! -e f.txt ]; then
   echo -e "#Population_1\tPopulation_2\tPopulation_3a\tPopualtion_3b\tf_hom\tf_d\tf" > f.txt
fi
if ! grep -q -P "concolor\troumanicus\teastern\twestern\teuropaeus" f.txt; then
   # In principle, P3a is the one expected to have contributed some variation to P2, and P3b is the one
   # that substitutes P2 to represent the situation of complete admixture. So, P3a shouls be the eastern
   # europaeus:
   R -q --no-save <estimate_f.R --args combined.tsv concolor roumanicus eastern western europaeus | grep -vP "^[#>\+]" >> f.txt
fi
if ! grep -q -P "concolor\troumanicus\twestern\teastern\teuropaeus" f.txt; then
   # And this is expected to underestimate the excess of ABBA between P2 and P3a, and therefore to
   # underestimate f.
   R -q --no-save <estimate_f.R --args combined.tsv concolor roumanicus western eastern europaeus | grep -vP "^[#>\+]" >> f.txt
fi

# Note that there are other estimators of the admixed portion of the genome, like f_dM. See:
#
#    https://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/
#
# and also the implementation in python by Simon Martin in his genomics.py module. We could just run its pipeline.
