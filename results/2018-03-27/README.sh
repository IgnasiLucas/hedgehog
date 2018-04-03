#!/bin/bash
#
#				2018-03-27
#				----------
#
# We received more sequence data, and Kristýna has already aligned it and
# merged it with the previous data. Here, I want to compare the completeness
# of the vcf files before and after the last sequencing run.
#
# The last vcf build before the last sequencing run is ../2017-05-25/hhog.vcf.
# It was obtained with freebayes. It includes a binary flag indicating what
# samples each allele is present in. It is not filtered by quality.

BEFORE=../2017-05-25/hhog.vcf

# The vcf obtained by Kristýna after the addition of the last batch of sequences
# is in /data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz. It was
# also obtained with freebayes, with similar options.

AFTER=/data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz
AFTERBAM=/data/kristyna/hedgehog/data_2018/merged_bam

# For the comparison to be accurate and fast, I could use only a subset of sites.
# For example, the sites covered in the only Hemichinus sample, which are the
# sites among which the abba/baba test is feasible.

if [ ! -e popmap.txt ]; then
   cp ../2017-05-25/populations.txt ./popmap.txt
fi

# The order of samples in the vcf files is not the same. I will select the sites
# where Hemichinus was covered in the vcf file before the last sequencing event,
# and the quality of the variant was at least 50.

HEMICHINUS_POSITION_BEFORE=$( head -n 6000 $BEFORE | gawk '(/^#CHROM/){for (i = 0; i <= NF - 10; i++) {if ($(i + 10) == "Er65_IS25") print i}}' )
echo "Hemichinus is sample number $HEMICHINUS_POSITION_BEFORE"
HEMICHINUS_FLAG=$( echo "2 ^ $HEMICHINUS_POSITION_BEFORE" | bc )
echo "The corresponding binary number is $HEMICHINUS_FLAG"

if [ ! -e comparison.txt ]; then
   if [ ! -e hemichinus_sites.txt ]; then
      gawk -f filtervcf.awk -v REQUIRED=$HEMICHINUS_FLAG $BEFORE | gawk '((/^[^#]/) && ($6 >= 50.0)){print $1 "\t" $2}' > hemichinus_sites.txt
   fi

   if [ ! -e before.recode.vcf ] && [ ! -e common_before.recode.vcf ]; then
      vcftools --vcf $BEFORE \
               --out before \
               --positions hemichinus_sites.txt \
               --recode 1> before.log 2> before.err &
   fi

   if [ ! -e after.recode.vcf ] && [ ! -e common_after.recode.vcf ]; then
      vcftools --gzvcf $AFTER \
               --out after \
               --positions hemichinus_sites.txt \
               --recode 1> after.log 2> after.err &
   fi

   wait

   # In order to use exactly the same sites, I will reduce the previous files
   # further. First, I need to exclude chromosomes present in only one file.
   # I manually determined them and use them below.

   if [ ! -e beforeafter.diff.sites_in_files ]; then
      vcftools --vcf before.recode.vcf \
               --out beforeafter \
               --not-chr NW_006805355.1 \
               --not-chr NW_006805505.1 \
               --not-chr NW_006805518.1 \
               --not-chr NW_006805635.1 \
               --not-chr NW_006805650.1 \
               --not-chr NW_006805660.1 \
               --not-chr NW_006805683.1 \
               --not-chr NW_006805778.1 \
               --not-chr NW_006805784.1 \
               --not-chr NW_006805791.1 \
               --not-chr NW_006805810.1 \
               --not-chr NW_006805815.1 \
               --not-chr NW_006805832.1 \
               --not-chr NW_006805936.1 \
               --not-chr NW_006806027.1 \
               --not-chr NW_006806028.1 \
               --not-chr NW_006806299.1 \
               --not-chr NW_006806372.1 \
               --not-chr NW_006806396.1 \
               --not-chr NW_006806515.1 \
               --not-chr NW_006806796.1 \
               --not-chr NW_006807000.1 \
               --not-chr NW_006807010.1 \
               --not-chr NW_006807104.1 \
               --not-chr NW_006807718.1 \
               --not-chr NW_006807727.1 \
               --not-chr NW_006807972.1 \
               --not-chr NW_006807985.1 \
               --not-chr NW_006808106.1 \
               --not-chr NW_006808186.1 \
               --not-chr NW_006808370.1 \
               --not-chr NW_006809027.1 \
               --not-chr NW_006809075.1 \
               --not-chr NW_006809618.1 \
               --diff after.recode.vcf \
               --diff-site
   fi

   if [ ! -e common_sites.txt ]; then
      grep B beforeafter.diff.sites_in_files | cut -f 1,2 > common_sites.txt
   fi

   if [ ! -e common_before.recode.vcf ]; then
      vcftools --vcf before.recode.vcf \
               --out common_before \
               --positions common_sites.txt \
               --recode 1> z1.log 2> z2.err &
      rm before.recode.vcf
   fi

   if [ ! -e common_after.recode.vcf ]; then
      vcftools --vcf after.recode.vcf \
               --out common_after \
               --positions common_sites.txt \
               --recode 1> z3.log 2> z4.err &
      rm after.recode.vcf
   fi

   wait

   # Now I want to create a file with at least six columns: chromosome,
   # site, variant quality before addition of data, variant quality
   # afterwards, number of samples with data before, and after data
   # addition.

   gawk '(/^[^#]/){
      WITHDATA = 0
      for (i = 10; i <= NF; i++) {
         if ($i ~ /[0-9]/) WITHDATA++
      }
      print $1 "\t" $2 "\t" $6 "\t" WITHDATA
   }' common_before.recode.vcf > z1 &
   gawk '(/^[^#]/){
      WITHDATA = 0
      for (i = 10; i <= NF; i++) {
         if ($i ~ /[0-9]/) WITHDATA++
      }
      print $6 "\t" WITHDATA
   }' common_after.recode.vcf > z2 &
   wait
   paste z1 z2 > comparison.txt
   rm z1 z2
fi

if [ -e common_before.recode.vcf ]; then rm common_before.recode.vcf; fi
if [ -e common_after.recode.vcf ];  then rm common_after.recode.vcf;  fi

# Visual inspection of the comparison file, with R, shows that there
# is substantial improvement. The average number of individuals with
# data per site changed from 32.6 to 40.5. The median changed from
# 37 to 44. The qualities also improved.
#
# However, some sites have fewer samples with data and/or lower qualities
# after the addition of sequencing data. A list of possible reasons for that:
#
#   * The filter settings in freebayes were different to produce the two original
#     vcf files. Before the addition of sequence data, the minimum mapping quality
#     was 30 and the minimum base quality, 20. After the addition of sequence
#     data, both thresholds were set to 15. Thus, part of the increase in the
#     number of samples with data may be due to the less stringent use of the
#     available reads. That could also explain some lower qualities. But it does
#     not help explain why in some sites there are fewer samples with data.
#
#   * The mapping of reads may have been different. I need to check that.
#
#
# On 2017-05-25, I plotted a histogram of the number of samples with data per site
# for the hhog.vcf file (before the additional sequencing). I want to plot the
# same now for both before and after.

if [ ! -e SamplesWithData.png ]; then
   if [ ! -e WithDataBefore.txt ]; then
      echo -e "#withData\tFrequency_before" > WithDataBefore.txt
      gawk '(/^[^#]/){
         split($8, INFO, /;/)
         for (i in INFO) {
            if (INFO[i] ~ /^NS=/) {
               F[substr(INFO[i],4)]++
            }
         }
      }END{
         for (f = 1; f <= 50; f++) print f "\t" F[f] + 0
      }' $BEFORE >> WithDataBefore.txt
   fi

   if [ ! -e WithDataAfter.txt ]; then
      echo "Frequency_after" > WithDataAfter.txt
      gunzip -c $AFTER | \
      gawk '(/^[^#]/){
         split($8, INFO, /;/)
         for (i in INFO) {
            if (INFO[i] ~ /^NS=/) {
               F[substr(INFO[i],4)]++
            }
         }
      }END{
         for (f = 1; f <= 50; f++) print F[f] + 0
      }' >> WithDataAfter.txt
   fi
   paste WithDataBefore.txt WithDataAfter.txt > SamplesWithData.txt
   rm WithDataBefore.txt WithDataAfter.txt
   gnuplot < SamplesWithData.gnp
fi

# The distribution of the number of samples with data per site has improved a lot.
# I need to know how many more reads were produced in the last sequencing run. An
# analysis similar to the one from 2016-11-30 is in order. There, I estimated the
# coverage per locus from the bam files, using bedtools. That way, the coverage
# values refer to loci, not to variable sites. The main difference between the
# distribution of number of samples with data per site and per locus is that the
# loci covered in only one or few individuals seem to include a disproportionate
# amount of variable sites. This is an indication of spurious mapping of some
# reads, probably due to low base qualities.
#
#                              First&Second	    Third	    Total
# Total number of read pairs	   ---  	241866527	   ---
# Demultiplexed             	319934047	220662575	540596622
# Mapped                        227639522	149124470	376763992
# Nuclear                   	

SAMPLE=(Er26_JUG4 Er27_SK32 Er28_112   Er29_122  Er30_79    Er31_453  Er32_183   Er33_211  Er34_197  Er35_209
        Er36_SK24 Er37_SK27 Er38_LB1   Er39_PL1  Er40_M2    Er41_GR36 Er42_GR35  Er43_SL7  Er44_VOJ1 Er45_BLG3
        Er46_RMN7 Er47_CR4  Er48_BH16  Er49_GR38 Er50_R3    Er51_436  Er52_451   Er53_ASR7 Er54_AU1  Er55_AU7
        Er56_AZ5  Er57_COR4 Er58_FI7   Er59_FR1  Er60_GR5   Er61_GR87 Er62_GR95  Er63_IR6  Er64_IS1  Er65_IS25
        Er66_IT3  Er67_IT5  Er68_PRT1B Er69_R2   Er70_RMN42 Er71_SAR2 Er72_SIE1B Er73_SNG1 Er74_SP16 Er75_TRC2A)

for i in `seq 0 24`; do
   if [ ! -e ${SAMPLE[$i]}.bed ]; then
      samtools view -f 64 -F 260 -b -u $AFTERBAM/${SAMPLE[$i]}'_msnc.bam' | bedtools merge -d 10 -c 3 -o count -i stdin > ${SAMPLE[$i]}.bed &
   fi
done
wait
for i in `seq 25 49`; do
   if [ ! -e ${SAMPLE[$i]}.bed ]; then
      samtools view -f 64 -F 260 -b -u $AFTERBAM/${SAMPLE[$i]}'_msnc.bam' | bedtools merge -d 10 -c 3 -o count -i stdin > ${SAMPLE[$i]}.bed &
   fi
done
wait

if [ ! -e pooled.bed ]; then
   samtools merge -u - $AFTERBAM/*_msnc.bam | samtools view -f 64 -F 260 -bu - | bedtools merge -d 10 -c 3 -o count -i stdin > pooled.bed
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

#if [ ! -e pooled_filtered.bed ]; then
#   gawk '(($4 >= 100) && ($4 <= 750)){print $1 "\t" $2 "\t" $3 "\t" $4}' pooled.bed > pooled_filtered.bed
#fi

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
# NW_006803924.1        9262    9457    511     Er54_AU1        NW_006803924.1  9310    9430    7       120
# NW_006803924.1        9262    9457    511     Er65_IS25       NW_006803924.1  9310    9457    281     147
# NW_006803924.1        10586   10908   208     Er54_AU1        NW_006803924.1  10672   10908   3       236
# NW_006803924.1        10586   10908   208     Er65_IS25       NW_006803924.1  10586   10896   27      310
# NW_006803924.1        11088   11333   460     Er54_AU1        NW_006803924.1  11212   11332   2       120
# NW_006803924.1        11088   11333   460     Er65_IS25       NW_006803924.1  11088   11332   266     244
# NW_006803924.1        11734   11992   325     Er54_AU1        NW_006803924.1  11734   11970   14      236
# NW_006803924.1        19108   19344   352     Er54_AU1        NW_006803924.1  19108   19344   12      236
# NW_006803924.1        19108   19344   352     Er65_IS25       NW_006803924.1  19108   19344   7       236
# NW_006803924.1        26318   26567   291     Er54_AU1        NW_006803924.1  26318   26554   9       236
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
# it below, counting on each site how many samples cover it with at least 1, 2... or
# 10 reads:
for i in 1 2 3 4 5 6 7 8 9 10 ; do
   if [ ! -e summary_coverage_min$i.txt ]; then
      gawk -v MIN=$i '(NR > 1){
         N=0
         for (i=4;i<=NF;i++) {
            if ($i >= (MIN + 0)) N++
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
      }'> summary_coverage_min$i.txt
   fi
done
