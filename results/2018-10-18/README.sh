#!/bin/bash
#
#				2018-10-18
#				==========
#
# Once I am familiar with the D-test (2018-10-15) and I know what parts of the
# genome of the admixed Er37_SK27 individual have an E. europaeus origin (2018-10-03),
# it's time to determine if the genomic regions where Er37_SK27 keeps E. europaeus
# ancestry are also characterized by a strong introgression signal among its peers.
#
# The idea may be too simple. Er37_SK27 must be the offspring of one normal and
# one hybrid, fully heterozygous parent. Hence, its approximately 25% of E. europaeus
# ancestry. Because Er37_SK27 survived development, and to the moment of sample
# collection, it can be said that its specific combination of ancestyries is not
# lethal, at least. And, because it survived that long, the most likely scenario is
# that the combination of ancestries is not deleterious. If natural selection removes
# from the population deleterious admixed genomic regions, then we expect signatures
# of admixture to be localized only where natural selections allows it. If that was
# most of the genome, we would not expect much coincidence among independently admixed
# individuals. But, if admixture was deleterious in most of the genome, even individuals
# with low levels of admixture who are not related, would share the admixed blocks.
#
# So, what I want to test is if the D statistic in the admixed regions is significantly
# higher than that outside the admixture regions. To assess significance, I need a null
# distribution... But let's check the difference first.

VCF=../2018-09-26/ErinaceusAndHemiechinus_r1e1c1h1.vcf
Er37=../2018-10-03/Er37_SK27/
FREQUENCIES=../2018-10-15/combined.tsv

if [ ! -e admixed_freqs.tsv ]; then
   if [ ! -e contig_lengths.txt ]; then
         gawk '(/^##contig=<ID=/){split($1,A,/[=,>]/); print A[3] "\t" A[5]}' $VCF > contig_lengths.txt
   fi
   # The files ancestry_lamp4.out, produced by LAMP, have as many rows as ancestries,
   # which are two. The first one is always 'E. roumanicus'. The second line is the one
   # that shows the amount of 'E. europaeus' ancestry for every SNP in the contig.
   # LAMP used 125299 sites, even though the Admixture analysis was run with fewer sites.
   # Admixture estimated 25.29% of europaeus ancestry. LAMP estimated 22.24% of SNPs or
   # 22.5 of all sites. Here, we use only 98650. The approach is to find the coordinates
   # of the ancestry blocks, and then impute ancestry to SNPs in a block that have not
   # been anlysed by either Admixture or LAMP. Actually, the proportion of sites in admixed
   # blocks is twice the proportion of europaeus ancestry (~47%), because all the europaeus
   # ancestry is heterozygous.
   if [ ! -e ancestry_blocks.bed ]; then
      touch ancestry_blocks.bed
      for contig in $(cut -f 1 contig_lengths.txt); do
         # Some contigs do not have results.
         if [ -s $Er37/$contig/gen000002/ancestry_lamp4.out ]; then
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
            }' $Er37/$contig/gen000002/ancestry_lamp4.out > z_$contig.txt
            # LAMP could have excluded some SNPs from analysis.
            if [ $(cat z_$contig.txt | wc -l) -eq $(cat $Er37/$contig/positions.txt | wc -l) ]; then
               paste $Er37/$contig/positions.txt z_$contig.txt | \
               gawk -v CONTIG=$contig -v LENGTH=$LENGTH 'BEGIN{
                  getline
                  START = 0
                  LAST = 1
                  TYPE = $2
               }($2 == TYPE){
                  LAST = $1
               }($2 != TYPE){
                  # I assign the recombination point to the midpoint between SNPs.
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
      gawk '(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' $FREQUENCIES > SNPs.bed
   fi
   bedtools intersect -wa -a SNPs.bed -b admixed.bed -sorted | cut -f 1,3 > admixed_SNP_positions.txt
   sed -r -i 's/$/\t/' admixed_SNP_positions.txt
   head -n 1 $FREQUENCIES > admixed_freqs.tsv
   grep -F -f admixed_SNP_positions.txt $FREQUENCIES >> admixed_freqs.tsv
#   rm admixed_SNP_positions.txt
fi

if [ ! -e roumanicus_freqs.tsv ]; then
   bedtools intersect -wa -a SNPs.bed -b roumanicus.bed -sorted | cut -f 1,3 > roumanicus_SNP_positions.txt
   sed -r -i 's/$/\t/' roumanicus_SNP_positions.txt
   head -n 1 $FREQUENCIES > roumanicus_freqs.tsv
   grep -F -f roumanicus_SNP_positions.txt $FREQUENCIES >> roumanicus_freqs.tsv
fi

if [ ! -e D.txt ]; then
   echo "#D statistic in regions where Er37_SK27 has mixed ancestries:" > D.txt
   R -q --no-save <abba_baba.R --args    admixed_freqs.tsv contig_lengths.txt concolor roumanicus europaeus 1e5 | grep -vP "^[#>\+]" >> D.txt
   echo "#D statistic in blocks where Er37_SK27 only has E. roumanicus ancestry:" >> D.txt
   R -q --no-save <abba_baba.R --args roumanicus_freqs.tsv contig_lengths.txt concolor roumanicus europaeus 1e5 | grep -vP "^[#>\+]" >> D.txt
fi

if [ ! -e f.txt ]; then
   echo "#f statistic in regions where Er37_SK27 has mixed ancestry" > f.txt
   R -q --no-save <estimate_f.R --args    admixed_freqs.tsv concolor roumanicus eastern western europaeus | grep -vP "^[#>\+]" >> f.txt
   echo "#f statistic in regions where E437_SK27 only has E. roumanicus ancestry:" >> f.txt
   R -q --no-save <estimate_f.R --args roumanicus_freqs.tsv concolor roumanicus eastern western europaeus | grep -vP "^[#>\+]" >> f.txt
fi

# Both D and f statistics are higher for regions where Er37_SK27 is admixed, as expected if natural
# selection is limiting the genomic regions where admixture is allowed. However, the difference is
# so small, that does not seem significant, according to the standard error of the D statistics.
# Bedtools offers a way to shuffle a genomic feature, in order to assess the significance of the
# association between features. Shuffling the admixed regions, and re-calculating D and t statistics
# would generate a null distribution of the D and t statistics, to test if the difference is significant.
# However, some of the admixed regions include almost a whole contig. Shuffling those will not produce
# enough variation unless I allow the shuffled block to extend beyond the end of a contig. That is
# not a good idea, because then the shuffled regions would not represent the same proportion of the
# genome. I think it is better to arbitrarily join contigs, to make fake scaffolds where shuffling
# is more effective. The easiest is to join them all together in a single scaffold, but there seems
# to be a limit in how large the chromosome can be. Let's try with two.

if [ ! -e FakeScaffold_lengths.txt ]; then
   cut -f 1 ancestry_blocks.bed | uniq > z_contigs_used.txt
   sed -r -i 's/$/\t/' z_contigs_used.txt
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

if [ ! -e Fake_admixed_freqs.tsv ]; then
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
   }' FakeScaffold_lengths.txt used_contig_lengths.txt admixed_freqs.tsv > Fake_admixed_freqs.tsv
fi

if [ ! -e Fake_roumanicus_freqs.tsv ]; then
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
   }' FakeScaffold_lengths.txt used_contig_lengths.txt roumanicus_freqs.tsv > Fake_roumanicus_freqs.tsv
fi

if [ ! -e FakeScaffold_freqs.tsv ]; then
   tail -n +2 Fake_admixed_freqs.tsv > z1
   tail -n +2 Fake_roumanicus_freqs.tsv > z2
   head -n  1 Fake_admixed_freqs.tsv > FakeScaffold_freqs.tsv
   cat z1 z2 | sort -nk 1,1 -k 2,2 >> FakeScaffold_freqs.tsv
fi

if [ ! -e Fake_SNPs.bed ]; then
   # This is just the coordinates of the SNPs with data, in bed format. I need it later
   # to intersect it with the shuffled blocks.
   gawk '(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' FakeScaffold_freqs.tsv > Fake_SNPs.bed
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

ORIGINAL=$(R -q --no-save <abba_baba_noJack.R --args      admixed_freqs.tsv concolor roumanicus europaeus | grep -vP "^[#>\+]")
    FAKE=$(R -q --no-save <abba_baba_noJack.R --args Fake_admixed_freqs.tsv concolor roumanicus europaeus | grep -vP "^[#>\+]")
echo "This two should be equal:"
echo $ORIGINAL
echo $FAKE

if [ ! -e null_diff_D.txt ] || [ ! -e null_diff_f.txt ]; then
   if [ ! -e null_D.txt ] || [ ! -e null_f.txt ]; then
      for i in $(seq 1 50); do
         ./get_nulls.sh $i &
      done
      wait

      for i in $(seq 1 50); do
         paste null_D1_$i.txt null_D2_$i.txt >> null_D.txt
         paste null_f1_$i.txt null_f2_$i.txt >> null_f.txt
         rm null_D1_$i.txt null_D2_$i.txt null_f1_$i.txt null_f2_$i.txt A_$i.* B_$i.* A_pos_$i.txt B_pos_$i.txt
      done
   fi

   gawk '{F[sprintf("%.3f", $1 - $2)]++}END{F["0.000"] += F["-0.000"]; delete F["-0.000"]; for (f in F) print f "\t" F[f]}' null_D.txt | sort -nk 1,1 > null_diff_D.txt
   echo -e "Diff\tf_hom\tf_d\tf" > null_diff_f.txt
   gawk '{
      F1[sprintf("%.3f", $6 - $14)]++
      F2[sprintf("%.3f", $7 - $15)]++
      F3[sprintf("%.3f", $8 - $16)]++
      Z[sprintf("%.3f", $6 - $14)]++
      Z[sprintf("%.3f", $7 - $15)]++
      Z[sprintf("%.3f", $8 - $16)]++
   }END{
      F1["0.000"] += F1["-0.000"]; delete F1["-0.000"]
      F2["0.000"] += F2["-0.000"]; delete F1["-0.000"]
      F3["0.000"] += F3["-0.000"]; delete F1["-0.000"]
      delete Z["-0.000"]
      for (z in Z) print z "\t" F1[z] + 0 "\t" F2[z] + 0 "\t" F3[z] + 0
   }' null_f.txt | sort -nk 1,1 >> null_diff_f.txt
fi

if [ ! -e summary.txt ]; then
   echo -e "# Values of abba/baba statistics in genomic regions with or without signals of"      > summary.txt
   echo -e "# recent admixture in Er37_SK27, which is excluded from the analysis. The p-values" >> summary.txt
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

# Values of abba/baba statistics in genomic regions with or without signals of
# recent admixture in Er37_SK27, which is excluded from the analysis. The p-values
# are computed from 5000 permutations of genomic regions.
#
# --------------------------------------------------------------------------------------
# Statistic	Admixed         	Not_admixed     	Difference	p-value
# -------------------------------------------------------------------------------------
# D       	0.159047740296386	0.133996097982671	0.0250516	0.2678
# f_hom   	0.0398418118213496	0.0349892634711391	0.00485255	0.3332
# f_d     	0.0235819111306138	0.0206065067310981	0.0029754	0.3432
# f     	0.0452036938062747	0.0410208465643631	0.00418285	0.3542
# --------------------------------------------------------------------------------------
