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
                  # We could assign the recombination point to the midpoint between SNPs...
#                 LAST = sprintf("%.0f", (LAST + $1) / 2)
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
   head -n 1 $FREQUENCIES > admixed_freqs.tsv
   grep -F -f admixed_SNP_positions.txt $FREQUENCIES >> admixed_freqs.tsv
#   rm admixed_SNP_positions.txt
fi

if [ ! -e roumanicus_freqs.tsv ]; then
   grep -v -F -f admixed_SNP_positions.txt $FREQUENCIES > roumanicus_freqs.tsv
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
   MAX_LENGTH=$(gawk '{S += $2}END{printf("%.0f\n", S/2)}' contig_lengths.txt)
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
   }' contig_lengths.txt > FakeScaffold_lengths.txt
fi

if [ ! -e FakeScaffold_freqs.tsv ]; then
   gawk 'BEGIN{
      FAKESCAF = 1
   }(FILENAME == "FakeScaffold_lengths.txt"){
      FAKELEN[$1] = $2
   }(FILENAME == "contig_lengths.txt"){
      if (CUMMUL == FAKELEN[FAKESCAF]){
         FAKESCAF++
         CUMMUL = 0
      }
      SCAF[$1] = FAKESCAF
      OFFSET[$1] = CUMMUL + 0
      CUMMUL += $2
   }(FILENAME ~ /.tsv$/){
      if (FNR == 1) print $0
      if (FNR > 1) {
         FAKEPOS = $2 + OFFSET[$1]
         $1 = SCAF[$1]
         $2 = FAKEPOS
         print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
      }
   }' FakeScaffold_lengths.txt contig_lengths.txt $FREQUENCIES > FakeScaffold_freqs.tsv
fi

if [ ! -e Fake_ancestry_blocks.bed ]; then
   gawk 'BEGIN{
      FAKESCAF = 1
   }(FILENAME == "FakeScaffold_lengths.txt"){
      FAKELEN[$1] = $2
   }(FILENAME == "contig_lengths.txt"){
      if (CUMMUL == FAKELEN[FAKESCAF]) {
         FAKESCAF++
         CUMMUL = 0
      }
      SCAF[$1] = FAKESCAF
      OFFSET[$1] = CUMMUL + 0
      CUMMUL += $2
   }(FILENAME == "ancestry_blocks.bed"){
      FAKESTART = $2 + OFFSET[$1]
      FAKEEND   = $3 + OFFSET[$1]
      $1 = SCAF[$1]
      $2 = FAKESTART
      $3 = FAKEEND
      print $1 "\t" $2 "\t" $3 "\t" $4
   }' FakeScaffold_lengths.txt contig_lengths.txt ancestry_blocks.bed > Fake_ancestry_blocks.bed
fi

if [ ! -e Fake_admixed.bed ]; then
   grep admixed Fake_ancestry_blocks.bed > Fake_admixed.bed
fi

if [ ! -e Fake_admixed_freqs.tsv ]; then
   if [ ! -e Fake_SNPs.bed ]; then
      gawk '(NR > 1){print $1 "\t" $2 - 1 "\t" $2}' FakeScaffold_freqs.tsv > Fake_SNPs.bed
   fi
   bedtools intersect -wa -a Fake_SNPs.bed -b Fake_admixed.bed -sorted | cut -f 1,3 > Fake_admixed_SNP_positions.txt
   head -n 1 FakeScaffold_freqs.tsv > Fake_admixed_freqs.tsv
   grep -F -f Fake_admixed_SNP_positions.txt FakeScaffold_freqs.tsv >> Fake_admixed_freqs.tsv
   rm Fake_admixed_SNP_positions.txt
fi

# I use a version of the abba_baba.R that does not run the jackknive method, which slows it down
# very much and is not necessary to generate a null distribution of D values.

ORIGINAL=$(R -q --no-save <abba_baba_noJack.R --args      admixed_freqs.tsv concolor roumanicus europaeus | grep -vP "^[#>\+]")
    FAKE=$(R -q --no-save <abba_baba_noJack.R --args Fake_admixed_freqs.tsv concolor roumanicus europaeus | grep -vP "^[#>\+]")
echo "This two should be equal:"
echo $ORIGINAL
echo $FAKE

if [ ! -e null_diff_D.txt ] || [ ! -e null_diff_f.txt ]; then
   for i in $(seq 1 10000); do
      # Shuffle the Fake_admixed.bed
      bedtools shuffle -noOverlapping -seed $i -i Fake_admixed.bed -g FakeScaffold_lengths.txt | sort -nk 1,1 -k 2,2 > A.bed
      # Get the its complement.
      bedtools complement -i A.bed -g FakeScaffold_lengths.txt > B.bed
      # Get the corresponding frequencies files.
      bedtools intersect -wa -a Fake_SNPs.bed -b A.bed -sorted | cut -f 1,3 > A_pos.txt
      bedtools intersect -wa -a Fake_SNPs.bed -b B.bed -sorted | cut -f 1,3 > B_pos.txt
      head -n 1 FakeScaffold_freqs.tsv | tee A.tsv > B.tsv
      grep -F -f A_pos.txt FakeScaffold_freqs.tsv >> A.tsv
      grep -F -f B_pos.txt FakeScaffold_freqs.tsv >> B.tsv
      # Estimate D and fs, and the differences.
      R -q --no-save <abba_baba_noJack.R --args A.tsv concolor roumanicus europaeus | grep -vP "^[#>\+]" >> null_D1.txt
      R -q --no-save <abba_baba_noJack.R --args B.tsv concolor roumanicus europaeus | grep -vP "^[#>\+]" >> null_D2.txt
      R -q --no-save <estimate_f.R --args A.tsv concolor roumanicus eastern western europaeus | grep -vP "^[#>\+]" >> null_f1.txt
      R -q --no-save <estimate_f.R --args B.tsv concolor roumanicus eastern western europaeus | grep -vP "^[#>\+]" >> null_f2.txt
   done
   paste null_D1.txt null_D2.txt | gawk '{F[sprintf("%.3f", $1 - $2)]++}END{for (f in F) print f "\t" F[f]}' | sort -nk 1,1 > null_diff_D.txt
   echo -e "Diff\tf_hom\tf_d\tf" > null_diff_f.txt
   paste null_f1.txt null_f2.txt | gawk '{
      F1[sprintf("%.3f", $6 - $14)]++
      F2[sprintf("%.3f", $7 - $15)]++
      F3[sprintf("%.3f", $8 - $16)]++
      Z[sprintf("%.3f", $6 - $14)]++
      Z[sprintf("%.3f", $7 - $15)]++
      Z[sprintf("%.3f", $8 - $16)]++
   }END{
      for (z in Z) print z "\t" F1[z] + 0 "\t" F2[z] + 0 "\t" F3[z] + 0
   }' | sort -nk 1,1 >> null_diff_f.txt
   rm A.bed B.bed A_pos.txt B_pos.txt A.tsv B.tsv
fi
