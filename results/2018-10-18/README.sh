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
   # blocks looks higher than 25%, presumably by chance.
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
   R -q --no-save <estimate_f.R --args   admixed_freqs.tsv concolor roumanicus western eastern europaeus | grep -vP "^[#>\+]" >> f.txt
   echo "#f statistic in regions where E437_SK27 only has E. roumanicus ancestry:" >> f.txt
   R -q --no-save <estimate_f.R --args roumanicus_freqs.tsv concolor roumanicus western eastern europaeus | grep -vP "^[#>\+]" >> f.txt
fi
