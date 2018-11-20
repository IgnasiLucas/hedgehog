#!/bin/bash
#
#            2018-10-24
#            ==========
#
# Er37_SK27, the hybrid individual, did not have its admixed genomic blocks associated
# with a significantly higher D or f statistic (2018-10-18). That was to be expected,
# since there was only one generation of selection since recombination of the F1. There
# is another admixed individual, Er55_AU7, which has inherited admixed blocks for some
# undetermined, longer number of generations. Natural selection had more time to remove
# deleterious introgressed haplotypes in its lineage, and it is worth checking if there
# is an association between the admixed blocks in Er55_AU7 and the D or t statistics.
#
# A difference with respect to the previous analysis is that the table of allele frequencies
# for the abba_baba analysis needs to be re-computed, because it should not include the
# individual that we are testing, Er55_AU7.
#
# Er55_AU7 has up to 2.37% of its genome admixed, according to LAMP and assuming 2500
# generations since admixture.

VCF=../2018-09-26/ErinaceusAndHemiechinus_r1e1c1h1.vcf
Er55=../2018-10-03/Er55_AU7/
ADMIXTURE_DIR=/data/kristyna/hedgehog/results_2018/05-06-2018/r10e10c4
FREQ2_AWK=../../bin/freq2.awk
POPDATA=../../data/populations.txt

if [ ! -e D.txt ] || [ ! -e f.txt ]; then
   if [ ! -e admixed_freqs.tsv ] || [ ! -e roumanicus_freqs.tsv ]; then
      echo $(date) " Calculating allele frequencies in admixed and non-admixed blocks."
      if [ ! -e frequencies.tsv ]; then
         # To calculate f statistics we need to partition the set of individuals in different ways,
         # with either the original P1, P2 and P3 populations, or in P1, P3a and P3b. This means that
         # a table of frequencies should include derived allele frequencies calculated for the two
         # alternative partitions. The script freq2.awk that I use to calculate allele frequencies
         # does not offer this option. Instead of updating it now, it is just convenient to run it
         # twice, with two different population maps: popmap1_test55.txt and popmap2_test55.txt.
         # then, the frequency tables will be combined for the common positions, just as in 2018-10-18.
         if [ ! -e popmap1_test55.txt ]; then
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
                  print $1 "\ttest"
               }
            }' $POPDATA | sort -k 2,2 > popmap1_test55.txt
         fi

         if [ ! -e popmap2_test55.txt ]; then
            if [ ! -e Q4.txt ]; then
               paste $ADMIXTURE_DIR/out.012.indv $ADMIXTURE_DIR/erinaceus_41_r10e10c4.4.Q > Q4.txt
            fi
            gawk 'BEGIN{
               POP[2] = "eastern"
               POP[3] = "concolor"
               POP[4] = "roumanicus"
               POP[5] = "western"
            }{
               POPULATION = "none"
               for (i = 2; i <= 5; i++) {
                  if ($i > 0.8) POPULATION = POP[i]
               }
               if (POPULATION == "none") POPULATION = "hybrid"
               if ($1 == "Er55_AU7") POPULATION = "test"
               print $1 "\t" POPULATION
            }END{
               print "Er36_SK24\troumanicus"
               print "Er58_FI7\teastern"      # This is a guess.
               print "Er61_GR87\troumanicus"
               print "Er62_GR95\troumanicus"
               print "Er67_IT5\teastern"      # This is a guess.
               print "Er65_IS25\tHemiechinus"
            }' Q4.txt | sort -k 2,2 > popmap2_test55.txt
         fi

         if [ ! -e creh1.tsv ]; then
            gawk -v P1="concolor" \
                 -v P2="roumanicus" \
                 -v P3="europaeus" \
                 -v OUTGROUP="Hemiechinus" \
                 -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MINOUT=1 \
                 -f $FREQ2_AWK popmap1_test55.txt $VCF | grep -v 999.9999 > creh1.tsv &
         fi

         if [ ! -e crweh1.tsv ]; then
            gawk -v P1="concolor" \
                 -v P2="roumanicus" \
                 -v P3="eastern" \
                 -v P4="western" \
                 -v OUTGROUP="Hemiechinus" \
                 -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MIN4=1 -v MINOUT=1 \
                 -f $FREQ2_AWK popmap2_test55.txt $VCF | grep -v 999.9999 > crweh1.tsv &
         fi
         wait
         cut -f 1,2 crweh1.tsv > zPositions.txt
         grep -F -f zPositions.txt creh1.tsv > z_creh1_subset.tsv
         paste crweh1.tsv z_creh1_subset.tsv | gawk '(($1 == $8) && ($2 == $9)){print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $12}' > frequencies.tsv
         rm zPositions.txt z_creh1_subset.tsv
      #  rm creh1.tsv crweh1.tsv popmap*
      fi

      if [ ! -e contig_lengths.txt ]; then
            gawk '(/^##contig=<ID=/){split($1,A,/[=,>]/); print A[3] "\t" A[5]}' $VCF > contig_lengths.txt
      fi
      # The files ancestry_lamp4.out, produced by LAMP, have as many rows as ancestries,
      # which are two. The first one is always 'E. roumanicus'. The second line is the one
      # that shows the amount of 'E. europaeus' ancestry for every SNP in the contig.
      # LAMP used 125299 sites, even though the Admixture analysis was run with fewer sites.
      if [ ! -e ancestry_blocks.bed ]; then
         touch ancestry_blocks.bed
         for contig in $(cut -f 1 contig_lengths.txt); do
            # Some contigs do not have results.
            if [ -s $Er55/$contig/gen002500/ancestry_lamp4.out ]; then
               LENGTH=$(grep $contig contig_lengths.txt | cut -f 2)
               gawk '(NR == 2){
                  for (i = 2; i <= NF; i++) {
                     if ($i == 0.0) {
                        print "roumanicus"
                     } else {
                        if (($i == 0.5) || ($i == 1.0)) {
                           print "admixed"
                        } else {
                           print $i
                        }
                     }
                  }
               }' $Er55/$contig/gen002500/ancestry_lamp4.out > z_$contig.txt
               # LAMP could have excluded some SNPs from analysis.
               if [ $(cat z_$contig.txt | wc -l) -eq $(cat $Er55/$contig/positions.txt | wc -l) ]; then
                  paste $Er55/$contig/positions.txt z_$contig.txt | \
                  gawk -v CONTIG=$contig -v LENGTH=$LENGTH 'BEGIN{
                     getline
                     START = 0
                     LAST = 1
                     TYPE = $2
                  }($2 == TYPE){
                     LAST = $1
                  }($2 != TYPE){
                     # Here I assign the recombination point to the midpoint between SNPs
                     LAST = sprintf("%.0f", (LAST + $1) / 2)
                     print CONTIG "\t" START "\t" LAST "\t" TYPE
                     START = LAST
                     LAST = START
                     TYPE = $2
                  }END{
                     LAST = LENGTH
                     print CONTIG "\t" START "\t" LAST "\t" TYPE
                  }' >> ancestry_blocks.bed
               else
                  echo "In contig $contig, the number of SNPs reported by LAMP is not the total number of SNPs."
               fi
               rm z_$contig.txt
            fi
         done
      fi

      if [ ! -e admixed.bed ]; then grep admixed ancestry_blocks.bed > admixed.bed; fi
      if [ ! -e roumanicus.bed ]; then grep roumanicus ancestry_blocks.bed > roumanicus.bed; fi
      if [ ! -e SNPs.bed ]; then
         gawk '(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' frequencies.tsv > SNPs.bed
      fi
      bedtools intersect -wa -a SNPs.bed -b admixed.bed -sorted | cut -f 1,3 > admixed_SNP_positions.txt
      # I need to add a tab at the end of each line to prevent matching undesired positions.
      sed -ri 's/$/\t/' admixed_SNP_positions.txt
      head -n 1 frequencies.tsv > admixed_freqs.tsv
      grep -F -f admixed_SNP_positions.txt frequencies.tsv >> admixed_freqs.tsv
   #  rm  admixed.bed roumanicus.bed
   fi

   if [ ! -e roumanicus_freqs.tsv ]; then
      # Note that the file frequencies.tsv contain more sites than those for which ancestry in Er55_AU7 is known.
      bedtools intersect -wa -a SNPs.bed -b roumanicus.bed -sorted | cut -f 1,3 > roumanicus_SNP_positions.txt
      sed -ri 's/$/\t/' roumanicus_SNP_positions.txt
      head -n 1 frequencies.tsv > roumanicus_freqs.tsv
      grep -F -f roumanicus_SNP_positions.txt frequencies.tsv >> roumanicus_freqs.tsv
   #  rm roumanicus_SNP_positions.txt
   fi

   if [ ! -e D.txt ]; then
      echo $(date) " Calculating D statistics."
      echo "#D statistic in regions where Er55_AU7 has mixed ancestries:" > D.txt
      R -q --no-save <abba_baba.R --args    admixed_freqs.tsv contig_lengths.txt concolor roumanicus europaeus 1e5 | grep -vP "^[#>\+]" >> D.txt
      echo "#D statistic in blocks where Er55_AU7 only has E. roumanicus ancestry:" >> D.txt
      R -q --no-save <abba_baba.R --args roumanicus_freqs.tsv contig_lengths.txt concolor roumanicus europaeus 1e5 | grep -vP "^[#>\+]" >> D.txt
   fi

   if [ ! -e f.txt ]; then
      echo $(date) " Calculating f statistics."
      echo "#f statistic in regions where Er55_AU7 has mixed ancestry" > f.txt
      R -q --no-save <estimate_f.R --args    admixed_freqs.tsv concolor roumanicus eastern western europaeus | grep -vP "^[#>\+]" >> f.txt
      echo "#f statistic in regions where Er55_AU7 only has E. roumanicus ancestry:" >> f.txt
      R -q --no-save <estimate_f.R --args roumanicus_freqs.tsv concolor roumanicus eastern western europaeus | grep -vP "^[#>\+]" >> f.txt
   fi
fi

if [ ! -e null_diff_D.txt ] || [ ! -e null_diff_f.txt ]; then
   echo $(date) " Generating the null distribution of differences in D and f between kinds of ancestry blocks."
   if [ ! -e FakeScaffold_lengths.txt ]; then
      cut -f 1 ancestry_blocks.bed | uniq > z_contigs_used.txt
      grep -F -f z_contigs_used.txt contig_lengths.txt > used_contig_lengths.txt
      rm z_contigs_used.txt
      MAX_LENGTH=$(gawk '{S += $2}END{printf("%.0f\n", S/2)}' used_contig_lengths.txt)
      gawk -v MAXLEN=$MAX_LENGTH 'BEGIN{
         SCAF = 1
      }{
         if (LEN < MAXLEN) {
            LEN += $2
         } else {
            print SCAF "\t" LEN
            LEN = $2
            SCAF++
         }
      }END{
         if (LEN <= MAXLEN) print SCAF "\t" LEN
      }' used_contig_lengths.txt > FakeScaffold_lengths.txt
   fi

   # I need to translate both the admixed_freqs.tsv file and also the whole set of allele frequencies
   # into the coordinates of the two fake scaffolds. However, the whole set should refer only to the
   # scaffolds actually used, namely the ones for which ancestry information of Er55_AU7 is available.
   # Thus, I can either translate admixed_freqs.tsv and roumanicus_freqs.tsv, and then join them, or
   # I could translate the subset of frequencies.tsv contained in used contigs, then translate the
   # ancestry_blocks.bed file, split it by ancestry, and intersect the result with the translated
   # frequencies.tsv to get translated versions of admixed_freqs.tsv and roumanicus_freqs.tsv. The
   # first option seems simpler. But note that I also need the translated version of ancestry_blocks.bed.

   if [ ! -e Fake_admixed.tsv ]; then
      gawk 'BEGIN{
         FAKE_SCAF = 1
      }(FILENAME == "FakeScaffold_lengths.txt"){
         FAKE_LEN[$1] = $2
      }(FILENAME == "used_contig_lengths.txt"){
         if (CUMMUL == FAKE_LEN[FAKE_SCAF]){
            FAKE_SCAF++
            CUMMUL = 0
         }
         SCAF[$1] = FAKE_SCAF
         OFFSET[$1] = CUMMUL + 0
         CUMMUL += $2
      }(FILENAME == "admixed_freqs.tsv"){
         if (FNR == 1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
         if (FNR > 1) {
            print SCAF[$1] "\t" $2 + OFFSET[$1] "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
         }
      }' FakeScaffold_lengths.txt used_contig_lengths.txt admixed_freqs.tsv > Fake_admixed.tsv
   fi
   if [ ! -e Fake_roumanicus.tsv ]; then
      gawk 'BEGIN{
         FAKE_SCAF = 1
      }(FILENAME == "FakeScaffold_lengths.txt"){
         FAKE_LEN[$1] = $2
      }(FILENAME == "used_contig_lengths.txt"){
         if (CUMMUL == FAKE_LEN[FAKE_SCAF]){
            FAKE_SCAF++
            CUMMUL = 0
         }
         SCAF[$1] = FAKE_SCAF
         OFFSET[$1] = CUMMUL + 0
         CUMMUL += $2
      }(FILENAME == "roumanicus_freqs.tsv"){
         if (FNR == 1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
         if (FNR > 1) {
            print SCAF[$1] "\t" $2 + OFFSET[$1] "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
         }
      }' FakeScaffold_lengths.txt used_contig_lengths.txt roumanicus_freqs.tsv > Fake_roumanicus.tsv
   fi
   if [ ! -e Fake_frequencies.tsv ]; then
      tail -n +2 Fake_admixed.tsv > z1
      tail -n +2 Fake_roumanicus.tsv > z2
      head -n  1 Fake_admixed.tsv > Fake_frequencies.tsv
      cat z1 z2 | sort -nk 1,1 -k 2,2 >> Fake_frequencies.tsv
   fi

   if [ ! -e Fake_SNPs.bed ]; then
      # This is just the coordinates of the SNPs with data, in bed format. I need it later,
      # in get_nulls.sh, to intersect it with the shuffled blocks.
      gawk '(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' Fake_frequencies.tsv > Fake_SNPs.bed
   fi

   if [ ! -e Fake_ancestry_blocks.bed ]; then
      gawk 'BEGIN{
         NEWSCAF = 1
         getline <"ancestry_blocks.bed"
         FAKESCAF = 1
         FAKESTART = 0
         FAKEEND = $3
         TYPE = $4
      }(FILENAME == "FakeScaffold_lengths.txt"){
         FAKELEN[$1] = $2
      }(FILENAME == "used_contig_lengths.txt"){
         if (CUMMUL == FAKELEN[NEWSCAF]) {
            NEWSCAF++
            CUMMUL = 0
         }
         SCAF[$1] = NEWSCAF
         OFFSET[$1] = CUMMUL + 0
         CUMMUL += $2
      }(FILENAME == "ancestry_blocks.bed"){
         if ((SCAF[$1] == FAKESCAF) && ($4 == TYPE)) {
            FAKEEND = $3 + OFFSET[$1]
         } else {
            print FAKESCAF "\t" FAKESTART "\t" FAKEEND "\t" TYPE
            FAKESCAF = SCAF[$1]
            FAKESTART = $2 + OFFSET[$1]
            FAKEEND = $3 + OFFSET[$1]
            TYPE = $4
         }
      }END{
         if ((SCAF[$1] == FAKESCAF) && ($4 == TYPE)) {
            print FAKESCAF "\t" FAKESTART "\t" FAKEEND "\t" TYPE
         }
      }' FakeScaffold_lengths.txt used_contig_lengths.txt ancestry_blocks.bed > Fake_ancestry_blocks.bed
   fi

   if [ ! -e Fake_admixed.bed ]; then
      grep admixed Fake_ancestry_blocks.bed > Fake_admixed.bed
   fi

   # I use a version of the abba_baba.R that does not run the jackknive method, which slows it down
   # very much and is not necessary to generate a null distribution of D values.

   ORIGINAL=$(R -q --no-save <abba_baba_noJack.R --args admixed_freqs.tsv concolor roumanicus europaeus | grep -vP "^[#>\+]")
       FAKE=$(R -q --no-save <abba_baba_noJack.R --args  Fake_admixed.tsv concolor roumanicus europaeus | grep -vP "^[#>\+]")
   echo "This two should be equal:"
   echo $ORIGINAL
   echo $FAKE
   if [ ! -e null_D.txt ] || [ ! -e null_f.txt ]; then
      for i in $(seq 1 50); do
         ./get_nulls.sh $i &
      done
      wait
      for i in $(seq 1 50); do
         paste null_D1_$i.txt null_D2_$i.txt >> null_D.txt
         paste null_f1_$i.txt null_f2_$i.txt >> null_f.txt
      done
      rm null_D1_*.txt null_D2_*.txt null_f1_*.txt null_f2_*.txt A_* B_*
   fi

   gawk '{F[sprintf("%.2f", $1 - $2)]++}END{F["0.00"] += F["-0.00"]; delete F["-0.00"]; for (f in F) print f "\t" F[f]}' null_D.txt | sort -nk 1,1 > null_diff_D.txt
   echo -e "Diff\tf_hom\tf_d\tf" > null_diff_f.txt
   gawk '{
      F1[sprintf("%.2f", $6 - $14)]++
      F2[sprintf("%.2f", $7 - $15)]++
      F3[sprintf("%.2f", $8 - $16)]++
      Z[sprintf("%.2f", $6 - $14)]++
      Z[sprintf("%.2f", $7 - $15)]++
      Z[sprintf("%.2f", $8 - $16)]++
   }END{
      F1["0.00"] += F1["-0.00"]; delete F["-0.00"]
      F2["0.00"] += F2["-0.00"]; delete F["-0.00"]
      F3["0.00"] += F3["-0.00"]; delete F["-0.00"]
      delete Z["-0.00"]
      for (z in Z) print z "\t" F1[z] + 0 "\t" F2[z] + 0 "\t" F3[z] + 0
   }' null_f.txt | sort -nk 1,1 >> null_diff_f.txt
#  rm admixed_freqs.tsv roumanicus_freqs.tsv Fake*
fi

if [ ! -e summary.txt ]; then
   echo -e "# Values of abba/baba statistics in genomic regions with or without signals of"      > summary.txt
   echo -e "# recent admixture in Er55_AU7, which is excluded from the analysis. The p-values"  >> summary.txt
   echo -e "# are one-tailed, computed from 5000 permutations of genomic regions."              >> summary.txt
   echo -e "#"                                                                                  >> summary.txt
   echo -e "#Statistic\tAdmixed\tNot_admixed\tDifference\tp-value"                              >> summary.txt
   gawk '((FILENAME == "D.txt") && (FNR == 2)){
      D1 = $4
   }((FILENAME == "D.txt") && (FNR == 4)){
      D2 = $4
   }((FILENAME == "f.txt") && (FNR == 2)){
      FHOM1 = $6; FD1 = $7; F1 = $8
   }((FILENAME == "f.txt") && (FNR == 4)){
      FHOM2 = $6; FD2 = $7; F2 = $8
   }(FILENAME == "null_diff_D.txt"){
      SUM_D += $2
      if ($1 >= (D1 - D2)) NUM_D += $2
   }((FILENAME == "null_diff_f.txt") && (FNR > 1)){
      SUM_FHOM += $2; SUM_FD += $3; SUM_F += $4
      if ($1 >= (FHOM1 - FHOM2)) NUM_FHOM += $2
      if ($1 >= (FD1 - FD2)) NUM_FD += $3
      if ($1 >= (F1 - F2)) NUM_F += $4
   }END{
      print "D\t"     D1    "\t" D2    "\t" D1 - D2       "\t" NUM_D / SUM_D
      print "f_hom\t" FHOM1 "\t" FHOM2 "\t" FHOM1 - FHOM2 "\t" NUM_FHOM / SUM_FHOM
      print "f_d\t"   FD1   "\t" FD2   "\t" FD1 - FD2     "\t" NUM_FD / SUM_FD
      print "f\t"     F1    "\t" F2    "\t" F1 - F2       "\t" NUM_F / SUM_F
   }' D.txt f.txt null_diff_D.txt null_diff_f.txt >> summary.txt
fi

