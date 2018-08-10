#!/bin/bash
#
#				2018-08-01
#				----------
#
# Our last results of the Admixture analysis, using a more permissive
# filtering of the available SNPs revealed an unexpected 4.3% contribution of
# concolor ancestry into an E. romanicus individual from Macedonia, Er26_JUG4.
# This is extremely unlikely to be true (Barbora), and I want to know how this
# signal was detected.
#
# In the 2018-06-30 folder, I see 3543 loci genotyped in Er26_JUG4 have log odds
# higher than 3 of having at least one allele contributed by the concolor ancestry.
# Note that the concolor ancestry is correctly identified by Admixture with K=3,
# because of SNPs fixed since divergence, even though only 5 individuals are used
# to estimate allele frequencies.

LASTDIR=`pwd | sed 's/2018-08-01/2018-06-30/'`
K=3

if [ ! -e frequencies.txt ]; then
   echo -e "#Contig\tPosition\tEuropaeus\tRoumanicus\tConcolor" > frequencies.txt
   paste $LASTDIR/positions.txt $LASTDIR/K3.P >> frequencies.txt
fi

if [ ! -e Er26.profile ]; then
   ln $LASTDIR/Er26_JUG4.out Er26.profile
fi

if [ ! -e Er26.Q ]; then
   paste $LASTDIR/SampleNames.txt $LASTDIR/K3.Q | grep Er26_ > Er26.Q
fi

for i in `seq 1 $K`; do
   if [ ! -e ancestry_$i.txt ]; then
      paste $LASTDIR/SampleNames.txt $LASTDIR/K3.Q | gawk -v ANCESTRY=$i '($(ANCESTRY + 1) > 0.5){print $1}' | grep -vP "Er37|Er26" > ancestry_$i.txt
   fi
done

if [ ! -e PopulationNames.txt ]; then
   gawk '(FILENAME ~ /popmap.txt/){POP[$1]=$NF}(FILENAME ~ /SampleNames.txt/){print POP[$1]}' $LASTDIR/popmap.txt $LASTDIR/SampleNames.txt > PopulationNames.txt
fi

if [ ! -e LOD3_ordered.012 ]; then
   if [ ! -e LOD3.012 ]; then
      cp $LASTDIR/SampleNames.txt LOD3.012
      gawk '($5 > 3){print $1 "\t" $2 "\t"}' Er26.profile > LOD3_positions.txt
      # File frequencies.txt has one added header line, just as the 012 has one header column.
      # command 'nl' adds line numbers.
      nl frequencies.txt | grep -Ff LOD3_positions.txt | cut -f 1 > z_ln
      for locus in `cat z_ln`; do
         cut -f $locus $LASTDIR/genotypes.012 > z1
         paste LOD3.012 z1 > z2
         mv z2 LOD3.012
      done
      rm z1 z_ln
   fi
   grep -f ancestry_3.txt LOD3.012  > LOD3_ordered.012
   grep "Er26" LOD3.012            >> LOD3_ordered.012
   grep -f ancestry_2.txt LOD3.012 >> LOD3_ordered.012
   grep "Er37" LOD3.012            >> LOD3_ordered.012
   grep -f ancestry_1.txt LOD3.012 >> LOD3_ordered.012
fi

if [ ! -e LOD3_P.txt ]; then
   grep -Ff LOD3_positions.txt frequencies.txt > LOD3_P.txt
fi

if [ ! -e LOD3_EmpiricalFreq.txt ]; then
   gawk 'BEGIN{
      POP["Er29_122"]   = "europaeus";  POP["Er30_79"]   = "europaeus";  POP["Er31_453"]  = "europaeus";  POP["Er33_211"]   = "europaeus"
      POP["Er34_197"]   = "europaeus";  POP["Er51_436"]  = "europaeus";  POP["Er52_451"]  = "europaeus";  POP["Er54_AU1"]   = "europaeus"
      POP["Er56_AZ5"]   = "europaeus";  POP["Er57_COR4"] = "europaeus";  POP["Er63_IR6"]  = "europaeus";  POP["Er71_SAR2"]  = "europaeus"
      POP["Er72_SIE1B"] = "europaeus";  POP["Er74_SP16"] = "europaeus";  POP["Er26_JUG4"] = "roumanicus"; POP["Er27_SK32"]  = "roumanicus"
      POP["Er28_112"]   = "roumanicus"; POP["Er32_183"]  = "roumanicus"; POP["Er35_209"]  = "roumanicus"; POP["Er37_SK27"]  = "roumanicus"
      POP["Er39_PL1"]   = "roumanicus"; POP["Er40_M2"]   = "roumanicus"; POP["Er41_GR36"] = "roumanicus"; POP["Er42_GR35"]  = "roumanicus"
      POP["Er43_SL7"]   = "roumanicus"; POP["Er44_VOJ1"] = "roumanicus"; POP["Er45_BLG3"] = "roumanicus"; POP["Er46_RMN7"]  = "roumanicus"
      POP["Er47_CR4"]   = "roumanicus"; POP["Er48_BH16"] = "roumanicus"; POP["Er50_R3"]   = "roumanicus"; POP["Er55_AU7"]   = "roumanicus"
      POP["Er60_GR5"]   = "roumanicus"; POP["Er66_IT3"]  = "roumanicus"; POP["Er69_R2"]   = "roumanicus"; POP["Er70_RMN42"] = "roumanicus"
      POP["Er38_LB1"]   = "concolor";   POP["Er49_GR38"] = "concolor";   POP["Er53_ASR7"] = "concolor";   POP["Er64_IS1"]   = "concolor"
      POP["Er75_TRC2A"] = "concolor"
   }{
      for (i = 2; i <= NF; i++) {
         if ($i >= 0) {
            ALTFREQ[POP[$1]][i-1] += $i
            WITHDATA[POP[$1]][i-1]++
         }
      }
      NUM_LOCI = NF - 1
   }END{
      for (locus = 1; locus <= NUM_LOCI; locus++) {
         outstring = ""
         GLOBALFREQ = (ALTFREQ["europaeus"][locus] + ALTFREQ["roumanicus"][locus] + ALTFREQ["concolor"][locus]) / ( 2 * (WITHDATA["europaeus"][locus] + WITHDATA["roumanicus"][locus] + WITHDATA["concolor"][locus]))
         if (WITHDATA["europaeus"][locus] > 0) {
            if (GLOBALFREQ >  0.5) outstring = outstring sprintf(  "%.3f",       ALTFREQ["europaeus"][locus]  / (2 * WITHDATA["europaeus"][locus]))
            if (GLOBALFREQ <= 0.5) outstring = outstring sprintf(  "%.3f", 1.0 - ALTFREQ["europaeus"][locus]  / (2 * WITHDATA["europaeus"][locus]))
         } else {
            outstring = outstring "NA"
         }
         if (WITHDATA["roumanicus"][locus] > 0) {
            if (GLOBALFREQ >  0.5) outstring = outstring sprintf("\t%.3f",       ALTFREQ["roumanicus"][locus] / (2 * WITHDATA["roumanicus"][locus]))
            if (GLOBALFREQ <= 0.5) outstring = outstring sprintf("\t%.3f", 1.0 - ALTFREQ["roumanicus"][locus] / (2 * WITHDATA["roumanicus"][locus]))
         } else {
            outstring = outstring "\tNA"
         }
         if (WITHDATA["concolor"][locus] > 0) {
            if (GLOBALFREQ >  0.5) outstring = outstring sprintf("\t%.3f",       ALTFREQ["concolor"][locus] / (2 * WITHDATA["concolor"][locus]))
            if (GLOBALFREQ <= 0.5) outstring = outstring sprintf("\t%.3f", 1.0 - ALTFREQ["concolor"][locus] / (2 * WITHDATA["concolor"][locus]))
         } else {
            outstring = outstring "\tNA"
         }
         print outstring
      }
   }' LOD3.012 > LOD3_EmpiricalFreq.txt
fi

if [ ! -e AFS.png ]; then
   if [ ! -e EurSpectrum.txt ]; then
      gawk '(/^[^#]/){
         if ($3 <  0.5) EUR[sprintf("%.2f", $3)]++
         if ($3 >= 0.5) EUR[sprintf("%.2f", 1-$3)]++
         N++
      }END{
         for (FREQ = 0; FREQ <= 1; FREQ += 0.01) {
            if (EUR[sprintf("%.2f", FREQ)] > 0) print FREQ "\t" EUR[sprintf("%.2f", FREQ)] / N
         }
      }' frequencies.txt > EurSpectrum.txt
   fi
   if [ ! -e RouSpectrum.txt ]; then
      gawk '(/^[^#]/){
         if ($4 <  0.5) ROU[sprintf("%.2f", $4)]++
         if ($4 >= 0.5) ROU[sprintf("%.2f", 1-$4)]++
         N++
      }END{
         for (FREQ = 0; FREQ <= 1; FREQ += 0.01) {
            if (ROU[sprintf("%.2f", FREQ)] > 0) print FREQ "\t" ROU[sprintf("%.2f", FREQ)] / N
         }
      }' frequencies.txt > RouSpectrum.txt
   fi
   if [ ! -e ConSpectrum.txt ]; then
      gawk '(/^[^#]/){
         if ($5 <  0.5) CON[sprintf("%.2f", $5)]++
         if ($5 >= 0.5) CON[sprintf("%.2f", 1-$5)]++
         N++
      }END{
         for (FREQ = 0; FREQ <= 1; FREQ += 0.01) {
            if (CON[sprintf("%.2f", FREQ)] > 0) print FREQ "\t" CON[sprintf("%.2f", FREQ)] / N
         }
      }' frequencies.txt > ConSpectrum.txt
   fi
   gnuplot plotSpectra.gnp
fi

# The comparison between the ML allele frequencies obtained by Admixture
# for each ancestry group and my own estimates (paste LOD3_P.txt LOD3_EmpiricalFreq.txt)
# showed that the frequencies agreed for the most part. However, in 12% of the positions
# I calculated the exact complement of what Admixture was calculating. I noticed this
# happened only when the observed global allele frequency was exactly 0.5. Then, I changed
# the equality sign in the awk commands above that calculate the observed frequencies, so
# that the alternative allele frequency is the major allele frequency only if it is > 0.5,
# instead of ">=". Now, only four loci (among the 1493 with LOD > 3) Admixture and I
# calculate the frequencies of different alleles. In any case, this should not affect the
# AncestryProfile.py script. It just looked as a bug here because I wanted to make sure
# that the frequencies estimated by admixture make sense. And they do. That's the main point.

if [ ! -e LOD3_comparison.txt ]; then
   gawk -f check.awk LOD3_ordered.012 LOD3_P.txt > LOD3_comparison.txt
fi

# About the genotype of the JUG4 individual at those sites with very likely concolor
# contribution, for the most part, it is heterozygous, when roumanicus and concolor have
# different alleles almost fixed. Thus, it does look real.
