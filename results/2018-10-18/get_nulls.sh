#!/bin/bash
for i in $(seq 1 100); do
   # Shuffle the Fake_admixed.bed
   bedtools shuffle -noOverlapping -seed $(( $1 * 1000 + $i )) -i Fake_admixed.bed -g FakeScaffold_lengths.txt | sort -nk 1,1 -k 2,2 > A_$1.bed
   # Get the its complement.
   bedtools complement -i A_$1.bed -g FakeScaffold_lengths.txt > B_$1.bed
   # Get the corresponding frequencies files.
   bedtools intersect -wa -a Fake_SNPs.bed -b A_$1.bed -sorted | cut -f 1,3 > A_pos_$1.txt
   sed -ri 's/$/\t/' A_pos_$1.txt
   bedtools intersect -wa -a Fake_SNPs.bed -b B_$1.bed -sorted | cut -f 1,3 > B_pos_$1.txt
   sed -ri 's/$/\t/' B_pos_$1.txt
   head -n 1 FakeScaffold_freqs.tsv | tee A_$1.tsv > B_$1.tsv
   grep -F -f A_pos_$1.txt FakeScaffold_freqs.tsv >> A_$1.tsv
   grep -F -f B_pos_$1.txt FakeScaffold_freqs.tsv >> B_$1.tsv
   # Estimate D and fs, and the differences.
   R -q --no-save <abba_baba_noJack.R --args A_$1.tsv concolor roumanicus europaeus | grep -vP "^[#>\+]" >> null_D1_$1.txt
   R -q --no-save <abba_baba_noJack.R --args B_$1.tsv concolor roumanicus europaeus | grep -vP "^[#>\+]" >> null_D2_$1.txt
   R -q --no-save <estimate_f.R --args A_$1.tsv concolor roumanicus eastern western europaeus | grep -vP "^[#>\+]" >> null_f1_$1.txt
   R -q --no-save <estimate_f.R --args B_$1.tsv concolor roumanicus eastern western europaeus | grep -vP "^[#>\+]" >> null_f2_$1.txt
done
