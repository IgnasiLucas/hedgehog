#!/bin/bash
#
#				2016-11-30
#				----------
#
# I want to count the number of loci covered in each sample, and to estimate the
# overlap among samples. I will use bedtools. I will use only first reads, which
# are the ones anchored at the restriction sites. Note that the sequence shown in
# sam format is in the same strand as the reference; so that sometimes it's been
# reverse-complemented. That is to say that the restriction site where the read 1
# starts can be either on 5-prime of reads that match the reference strand, or on
# 3-prime of reads that originally matched the reverse complementary of the reference.
#

DATADIR=../../data
SAMPLE=(Er26_JUG4 Er27_SK32 Er28_112   Er29_122  Er30_79    Er31_453  Er32_183   Er33_211  Er34_197  Er35_209
        Er36_SK24 Er37_SK27 Er38_LB1   Er39_PL1  Er40_M2    Er41_GR36 Er42_GR35  Er43_SL7  Er44_VOJ1 Er45_BLG3
        Er46_RMN7 Er47_CR4  Er48_BH16  Er49_GR38 Er50_R3    Er51_436  Er52_451   Er53_ASR7 Er54_AU1  Er55_AU7
        Er56_AZ5  Er57_COR4 Er58_FI7   Er59_FR1  Er60_GR5   Er61_GR87 Er62_GR95  Er63_IR6  Er64_IS1  Er65_IS25
        Er66_IT3  Er67_IT5  Er68_PRT1B Er69_R2   Er70_RMN42 Er71_SAR2 Er72_SIE1B Er73_SNG1 Er74_SP16 Er75_TRC2A)

for i in `seq 0 24`; do
   if [ ! -e ${SAMPLE[$i]}.bed ]; then
      samtools view -f 64 -F 260 -b -u $DATADIR/${SAMPLE[$i]}'_sortednuc.bam' | bedtools merge -d 10 -c 3 -o count -i stdin > ${SAMPLE[$i]}.bed &
   fi
done
wait
for i in `seq 25 49`; do
   if [ ! -e ${SAMPLE[$i]}.bed ]; then
      samtools view -f 64 -F 260 -b -u $DATADIR/${SAMPLE[$i]}'_sortednuc.bam' | bedtools merge -d 10 -c 3 -o count -i stdin > ${SAMPLE[$i]}.bed &
   fi
done
wait

if [ ! -e pooled.bed ]; then
   samtools merge -u - $DATADIR/*sortednuc.bam | samtools view -f 64 -F 260 -bu - | bedtools merge -d 10 -c 3 -o count -i stdin > pooled.bed
fi

# Once I have the set of loci ever covered by a sample, I notice the distribution
# of coverage quite irregular. I want to plot it.
if [ ! -e coverage.png ]; then
   if [ ! -e coverage.txt ]; then
      gawk '{F[$4]++}END{for (f in F) print f "\t" F[f]}' pooled.bed | sort -nk 1 | \
      gawk 'BEGIN{print "Coverage\tFrequency\tReads"}{S += $1 * $2; print $1 "\t" $2 "\t" S}' > coverage.txt
   fi
   R --no-save < plot_coverage.R
fi

# Although some sites are covered too much, the proportion of reads lost in those
# repetitive loci is not worrisome. From the plot, it is clear that we can focus
# on sites covered overall between 100 and 750 times:
if [ ! -e pooled_filtered.bed ]; then
   gawk '(($4 >= 100) && ($4 <= 750)){print $1 "\t" $2 "\t" $3 "\t" $4}' pooled.bed > pooled_filtered.bed
fi

# I tried unionbedg, which generates a matrix comparing coverage among the files.
# However, it splits each record as many times as necessary to distinguish the pieces
# covered by some but not other samples. In a way, it spoils the merging of features
# done before and generates a too big and messy dataset. I would like somethin more
# compact, just telling for each of the intervals in pooled_filtered.bed if it is
# covered by a sample at all or not. I should use the intersect function.
#
# The intersect function will also generate several lines per interval, one for each
# sample overlapping the interval, like this:
#
# bedtools intersect -a pooled_filtered.bg -b Er54_AU1.bed Er62_GR95.bed Er65_IS25.bed -names Er54_AU1 Er62_GR95 Er65_IS25-sorted -wo | head
# NW_006803924.1	9262	9457	511	Er54_AU1	NW_006803924.1	9310	9430	7	120
# NW_006803924.1	9262	9457	511	Er65_IS25	NW_006803924.1	9310	9457	281	147
# NW_006803924.1	10586	10908	208	Er54_AU1	NW_006803924.1	10672	10908	3	236
# NW_006803924.1	10586	10908	208	Er65_IS25	NW_006803924.1	10586	10896	27	310
# NW_006803924.1	11088	11333	460	Er54_AU1	NW_006803924.1	11212	11332	2	120
# NW_006803924.1	11088	11333	460	Er65_IS25	NW_006803924.1	11088	11332	266	244
# NW_006803924.1	11734	11992	325	Er54_AU1	NW_006803924.1	11734	11970	14	236
# NW_006803924.1	19108	19344	352	Er54_AU1	NW_006803924.1	19108	19344	12	236
# NW_006803924.1	19108	19344	352	Er65_IS25	NW_006803924.1	19108	19344	7	236
# NW_006803924.1	26318	26567	291	Er54_AU1	NW_006803924.1	26318	26554	9	236
#
# Where the last column is the number of bases overlapping the interval. The gawk script
# below turns this list in a matrix. In order to take into account the number of bases that
# overlap, I will re-calculate each sample's depth of covarge at an interval, multiplying the
# number of reads that had been merged by the proportion of the length originally spanned by
# those reads that actually overlap: $9 * $10 / ($8 - $7).

if [ ! -e coverage_matrix.txt ]; then
   HEADER="#CHR\tSTART\tEND"`printf "\t%s" "${SAMPLE[@]}"`
   bedtools intersect -a pooled_filtered.bed -b `printf "%s.bed " "${SAMPLE[@]}"` -names "${SAMPLE[@]}" -sorted -wo | \
   gawk -v HEADER="$HEADER" 'function printline(POS, COV){
         print POS "\t" COV["Er26_JUG4"]  + 0 "\t" COV["Er27_SK32"] + 0 "\t" COV["Er28_112"]   + 0 "\t" COV["Er29_122"]  + 0 "\t" COV["Er30_79"]    + 0 "\t" COV["Er31_453"]   + 0 "\t" \
                        COV["Er32_183"]   + 0 "\t" COV["Er33_211"]  + 0 "\t" COV["Er34_197"]   + 0 "\t" COV["Er35_209"]  + 0 "\t" COV["Er36_SK24"]  + 0 "\t" COV["Er37_SK272"] + 0 "\t" \
                        COV["Er38_LB1"]   + 0 "\t" COV["Er39_PL1"]  + 0 "\t" COV["Er40_M2"]    + 0 "\t" COV["Er41_GR36"] + 0 "\t" COV["Er42_GR35"]  + 0 "\t" COV["Er43_SL7"]   + 0 "\t" \
                        COV["Er44_VOJ1"]  + 0 "\t" COV["Er45_BLG3"] + 0 "\t" COV["Er46_RMN7"]  + 0 "\t" COV["Er47_CR4"]  + 0 "\t" COV["Er48_BH16"]  + 0 "\t" COV["Er49_GR38"]  + 0 "\t" \
                        COV["Er50_R3"]    + 0 "\t" COV["Er51_436"]  + 0 "\t" COV["Er52_451"]   + 0 "\t" COV["Er53_ASR7"] + 0 "\t" COV["Er54_AU1"]   + 0 "\t" COV["Er55_AU7"]   + 0 "\t" \
                        COV["Er56_AZ5"]   + 0 "\t" COV["Er57_COR4"] + 0 "\t" COV["Er58_FI7"]   + 0 "\t" COV["Er59_FR1"]  + 0 "\t" COV["Er60_GR5"]   + 0 "\t" COV["Er61_GR87"]  + 0 "\t" \
                        COV["Er62_GR95"]  + 0 "\t" COV["Er63_IR6"]  + 0 "\t" COV["Er64_IS1"]   + 0 "\t" COV["Er65_IS25"] + 0 "\t" COV["Er66_IT3"]   + 0 "\t" COV["Er67_IT5"]   + 0 "\t" \
                        COV["Er68_PRT1B"] + 0 "\t" COV["Er69_R2"]   + 0 "\t" COV["Er70_RMN42"] + 0 "\t" COV["Er71_SAR2"] + 0 "\t" COV["Er72_SIE1B"] + 0 "\t" COV["Er73_SNG1"]  + 0 "\t" \
                        COV["Er74_SP16"]  + 0 "\t" COV["Er75_TRC2A"] + 0
   }BEGIN{
      print HEADER
   }(NR == 1){
      POS = $1 "\t" $2 "\t" $3
   }($1 "\t" $2 "\t" $3 != POS){
      printline(POS,COV)
      delete COV
      POS = $1 "\t" $2 "\t" $3
      COV[$5] = $9 * $10 / ($8 - $7)
   }($1 "\t" $2 "\t" $3 == POS){
      COV[$5] = $9 * $10 / ($8 - $7)
   }END{
      printline(POS,COV)
   }' > coverage_matrix.txt
fi

# The coverage_matrix.txt file is 95MB, which is kind of too big. I will summarize
# it below, counting on each site how many samples cover it with at least 3 reads:
if [ ! -e summary_coverage.txt ]; then
   gawk '(NR > 1){
      N=0
      for (i=4;i<=NF;i++) {
         if ($i >= 3) N++
      }
      F[N]++
   }END{
      for (f in F) print f "\t" F[f]
   }' coverage_matrix.txt | sort -nrk 1 | gawk 'BEGIN{
      print "#Num.Samples\tNum.Sites\tAccumulated"
      ACCUMULATED = 0
   }{
      ACCUMULATED += $2
      print $1 "\t" $2 "\t" ACCUMULATED
   }'> summary_coverage.txt
fi

# CONCLUSIONS
# -----------
#
# I did not count the original number of reads directly, but from the pooled.bed
# file, I estimate a total number of (mapped) reads in the order 273 millions, of
# which, around 248 million are used to cover sites that get a total depth between
# 100 and 750. There are 662655 such sites. Thus, if the reads were evenly distributed
# among samples (50) and sites, each site would get 7.5 reads. In reality, though,
# there are only 157661 sites covered by at least 40 samples. This is much less than
# the original number of sites, and the reduction may be explained by the variance
# in coverage per site. However, this seems much more than the number of variable
# sites covered by all individuals of either population 1 () or population 2 (). Now,
# the question is how do the loci covered with some reads compare to the variable
# sites identified in 2016-11-24. In particular, I would like to see the distribution
# of the number of variable sites per locus.
