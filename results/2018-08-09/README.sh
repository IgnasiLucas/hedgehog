#!/bin/bash
#
#				2018-08-09
#				==========
#
# The fact is that we did not filter potential PCR duplicates from our dataset.
# I want to check what impact we expect from them. First, I will quantify them.
# I read this blog entry:
#
#   https://www.molecularecologist.com/2016/08/the-trouble-with-pcr-duplicates/
#
# It explains that PCR duplicates in RADseq are always present, because PCR amplification
# is one of the steps in library preparation. If sequencing was paired-end, the
# mapped position of the second read can be used to identify duplicates, because
# sonication should produce different second reads in each template.
#
# The number or proportion of reads that are redundant is one concern. But I want
# to know how duplicates are distributed: if many are duplicates of the same few
# reads, or if every read may get a few duplicates.
#
# I use the bam files:

BAM2016=/data/kristyna/hedgehog/data_2016/bam
BAM2018=/data/kristyna/hedgehog/data_2018/bam
BAMROOT=/data/kristyna/hedgehog/data_
VCFFILE=../2018-07-25/flagged.vcf

SAMPLE=(Er26_JUG4 Er27_SK32 Er28_112   Er29_122  Er30_79    Er31_453  Er32_183   Er33_211  Er34_197  Er35_209
        Er36_SK24 Er37_SK27 Er38_LB1   Er39_PL1  Er40_M2    Er41_GR36 Er42_GR35  Er43_SL7  Er44_VOJ1 Er45_BLG3
        Er46_RMN7 Er47_CR4  Er48_BH16  Er49_GR38 Er50_R3    Er51_436  Er52_451   Er53_ASR7 Er54_AU1  Er55_AU7
        Er56_AZ5  Er57_COR4 Er58_FI7   Er59_FR1  Er60_GR5   Er61_GR87 Er62_GR95  Er63_IR6  Er64_IS1  Er65_IS25
        Er66_IT3  Er67_IT5  Er68_PRT1B Er69_R2   Er70_RMN42 Er71_SAR2 Er72_SIE1B Er73_SNG1 Er74_SP16 Er75_TRC2A)

# I select second reads, order them by position, and count the number of reads
# mapped to the same position, and make a histogram of those counts for each
# library.

for year in 2016 2018; do
   for sample in ${SAMPLE[@]}; do
      if [ ! -e $sample'_'$year'_dups.hist' ]; then
         samtools view -h -u -F 260 -f 128 $BAMROOT$year/bam/$sample'_'$year.bam | \
         samtools sort -O SAM -T $sample - | \
         gawk 'BEGIN{
            CONTIG = "NW_006803924.1"
            POSITION = 1
            N = 0
            MAX_N = 1
         }(/^[^@]/){
            if (($3 == CONTIG) && ($4 == POSITION)) {
               N++
            } else {
               if (N > MAX_N) MAX_N = N
               F[N]++
               N = 1
               CONTIG = $3
               POSITION = $4
            }
         }END{
            if (N > MAX_N) MAX_N = N
            F[N]++
            SUM = 0
            for (n = 1; n <= MAX_N; n++) {
               SUM += n * F[n]
               print n "\t" F[n] "\t" SUM
            }
         }' > $sample'_'$year'_dups.hist' && echo "#end" >> $sample'_'$year'_dups.hist'
         sleep 3
         # samtools seems to not finish some jobs. I add the sleep statement and the
         # echo to make sure it works.
      fi
   done
   wait
done

# The main concern about PCR duplicates pointed out in the blog mentioned above
# is that they may inflate the number of homozygous genotypes. In this sense, the
# effect of PCR duplicates is the same as that of low coverage. In bi-allelic variants
# this should not bias allele frequency estimates (see 2018-03-27b). Only inference
# based on genotype (rather than allele) frequencies would be affected. Neither the
# Admixture analysis, nor the RADpainter analysis should be affected, I think.
# In any case, the overall failure to call heterozygous genotypes can be assessed
# looking at the proportions of the two alternative alleles in genotypes with two
# alleles. If PCR bias is a general concern, one of the two allelels will generally
# appear in more than 50% of the read counts. This can be seen in the final vcf file.
# Actually, the information field AB is precisely the allele balance at heterozygous
# genotypes, and it is reported for each site. AB may be zero if no heterozygous
# individuals are seen at a variable site. Excluding those, I see the following
# distribution of AB in the 2018-07-25/flagged.vcf file:

if [ ! -e AB.hist ]; then
   gawk '(/^[^#]/){
      split($8,INFO,/;/)
      for (FIELD in INFO) {
         split(INFO[FIELD], TAGVAL, /=/)
         if (TAGVAL[1] == "AB") {
            if (TAGVAL[2] > 0) F[sprintf("%.2f", TAGVAL[2])]++
            break
         }
      }
   } END {
      for (f in F) print f "\t" F[f]
   }' $VCFFILE | sort -gk 1,1 > AB.hist
fi

# Heterozygous sites have on average 48% reference reads, which is close to the expectation,
# but there is actually considerable dispersion of the the allele balance.
# This suggests that PCR duplicates do affect the observed allele balance at heterozygous
# sites, to the point of masking some heterozygous genotypes.
#
# Because PCR duplicates are library-specific, it may be interesting to look at the allele
# balance per individual, across sites. The rationale here is that some individuals may
# have more unbalanced allele observations than others, while all individuals from the
# same population are expected to have similar genome-wide heterozygosities under random
# mating. Note, however, that all individuals have reads from different libraries, the
# first ones with high levels of duplicates, and subsequent ones with less (see below).
#
# The random deviation from the 50% expectation due to low coverage makes one of the alleles
# look more frequent then the other. The same happens with PCR duplicates. Because in both
# cases the allele that gets sampled more often at each site is randomly chosen, I should
# not see an overall deviation from the 50% balance between reference and alternative
# allele. Instead, I must identify the alleles by how frequently they are observed, as
# 'minor' and 'major', at each site. If only low coverage was causing dispersion, 50%
# would still be the average. PCR duplicates, however, will cause the average to drop.
# I can estimate the overall average frequency of the least observed allele by adding
# up allele observations across heterozygous sites. That will give an accurate mean.
# I can also limit the effect of low coverage by setting a minimum coverage.

MINCOV=5
if [ ! -e AB_perSample.hist ]; then
      gawk -v MINCOV=$MINCOV 'BEGIN{
         # This order was manually determined from the map to present results by species,
         # and more or less from east to west, within a species. It may not be accurate.
         ORDER[1]  = 40; ORDER[2]  = 36; ORDER[3]  = 21; ORDER[4]  = 23; ORDER[5]  = 13; ORDER[6]  = 41; ORDER[7]  = 39; ORDER[8]  = 26; ORDER[9]  = 33
         ORDER[10] = 32; ORDER[11] = 19; ORDER[12] = 38; ORDER[13] = 27; ORDER[14] = 18; ORDER[15] = 34; ORDER[16] = 44; ORDER[17] = 10; ORDER[18] = 28
         ORDER[19] = 12; ORDER[20] = 42; ORDER[21] = 25; ORDER[22] = 47; ORDER[23] = 46; ORDER[24] = 16; ORDER[25] = 49; ORDER[26] = 22; ORDER[27] = 24
         ORDER[28] = 43; ORDER[29] = 20; ORDER[30] = 17; ORDER[31] = 29; ORDER[32] = 30; ORDER[33] = 37; ORDER[34] = 48; ORDER[35] = 15; ORDER[36] = 50
         ORDER[37] = 35; ORDER[38] = 54; ORDER[39] = 11; ORDER[40] = 53; ORDER[41] = 31; ORDER[42] = 52; ORDER[43] = 45; ORDER[44] = 51; ORDER[45] = 14
         POPULATION["Er56_AZ5"]   = "europaeus";  POPULATION["Er63_IR6"]  = "europaeus";  POPULATION["Er30_79"]   = "europaeus";  POPULATION["Er74_SP16"]  = "europaeus";  POPULATION["Er52_451"]  = "europaeus"
         POPULATION["Er71_SAR2"]  = "europaeus";  POPULATION["Er57_COR4"] = "europaeus";  POPULATION["Er67_IT5"]  = "europaeus";  POPULATION["Er31_453"]   = "europaeus";  POPULATION["Er51_436"]  = "europaeus"
         POPULATION["Er72_SIE1B"] = "europaeus";  POPULATION["Er34_197"]  = "europaeus";  POPULATION["Er29_122"]  = "europaeus";  POPULATION["Er54_AU1"]   = "europaeus";  POPULATION["Er33_211"]  = "europaeus"
         POPULATION["Er58_FI7"]   = "europaeus";  POPULATION["Er66_IT3"]  = "roumanicus"; POPULATION["Er32_183"]  = "roumanicus"; POPULATION["Er43_SL7"]   = "roumanicus"; POPULATION["Er55_AU7"]  = "roumanicus"
         POPULATION["Er27_SK32"]  = "roumanicus"; POPULATION["Er48_BH16"] = "roumanicus"; POPULATION["Er47_CR4"]  = "roumanicus"; POPULATION["Er35_209"]   = "roumanicus"; POPULATION["Er60_GR5"]  = "roumanicus"
         POPULATION["Er26_JUG4"]  = "roumanicus"; POPULATION["Er40_M2"]   = "roumanicus"; POPULATION["Er44_VOJ1"] = "roumanicus"; POPULATION["Er41_GR36"]  = "roumanicus"; POPULATION["Er42_GR35"] = "roumanicus"
         POPULATION["Er39_PL1"]   = "roumanicus"; POPULATION["Er45_BLG3"] = "roumanicus"; POPULATION["Er28_112"]  = "roumanicus"; POPULATION["Er70_RMN42"] = "roumanicus"; POPULATION["Er46_RMN7"] = "roumanicus"
         POPULATION["Er69_R2"]    = "roumanicus"; POPULATION["Er50_R3"]   = "roumanicus"; POPULATION["Er36_SK24"] = "roumanicus"; POPULATION["Er53_ASR7"]  = "roumanicus"; POPULATION["Er37_SK27"] = "hybrid"
         POPULATION["Er49_GR38"]  = "concolor";   POPULATION["Er64_IS1"]  = "concolor";   POPULATION["Er38_LB1"]  = "concolor";   POPULATION["Er75_TRC2A"] = "concolor";   POPULATION["Er65_IS25"] = "hemiechinus"
      }(/^#CHROM/){
         for (i = 10; i <= NF; i++) SAMPLE[i] = $i
         NUMSAMPLES = NF - 9
      }(/^[^#]/){
         split($9,FORMAT,/:/)
         for (FIELD in FORMAT) {
            if (FORMAT[FIELD] == "AD") AD = FIELD
            if (FORMAT[FIELD] == "GT") GT = FIELD
         }
         for (i = 10; i <= NF; i++) {
            split($i, SAMPLEINFO, /:/)
            split(SAMPLEINFO[AD], ALLELECOUNT, /,/)
            COV = ALLELECOUNT[1] + ALLELECOUNT[2]
            if ((ALLELECOUNT[1] != ".") && (COV >= MINCOV)) {
               NUMSITES[SAMPLE[i]]++
               GTFREQ[SAMPLE[i]][SAMPLEINFO[GT]]++
               if (ALLELECOUNT[1] <= ALLELECOUNT[2]) {MINOR=1; MAJOR=2} else {MINOR=2; MAJOR=1}
               AB[SAMPLE[i]][sprintf("%.1f", ALLELECOUNT[MINOR] / COV)]++
               if (SAMPLEINFO[GT] == "0/1") {
                  MINORCOV[SAMPLE[i]] += ALLELECOUNT[MINOR]
                  TOTALCOV[SAMPLE[i]] += COV
               }
            }
         }
      } END {
         LINE = "#Population     "
         for (i = 1; i <= NUMSAMPLES; i++) {
            LINE = LINE "\t" POPULATION[SAMPLE[ORDER[i]]]
         }
         print LINE
         LINE = "#Sample         "
         for (i = 1; i <= NUMSAMPLES; i++) {
            LINE = LINE "\t" SAMPLE[ORDER[i]]
         }
         print LINE
         LINE = "#Heterozygosity1"
         for (i = 1; i <= NUMSAMPLES; i++) {
            LINE = LINE sprintf("\t%.4f", GTFREQ[SAMPLE[ORDER[i]]]["0/1"] / NUMSITES[SAMPLE[ORDER[i]]])
            # This is the heterozygosity as actually reported in the vcf file, counting the believed heterozygous sites.
         }
         print LINE
         LINE = "#Heterozygosity2"
         for (i = 1; i <= NUMSAMPLES; i++) {
            LINE = LINE sprintf("\t%.4f", (NUMSITES[SAMPLE[ORDER[i]]] - AB[SAMPLE[ORDER[i]]]["0.0"]) / NUMSITES[SAMPLE[ORDER[i]]])
            # This heterozygosity counts as heterozygous sites with at least 1 different allele observation in 10 reads.
         }
         print LINE
         LINE = "#Average_AB"
         for (i = 1; i <= NUMSAMPLES; i++) {
            if (TOTALCOV[SAMPLE[ORDER[i]]] > 0) {
               LINE = LINE sprintf("\t%.6f", MINORCOV[SAMPLE[ORDER[i]]] / TOTALCOV[SAMPLE[ORDER[i]]])
            } else {
               LINE = LINE "\tnan"
            }
         }
         print LINE
         for (ab = 0.0; ab <= 0.51; ab += 0.1) {
            LINE = sprintf("%.1f", ab)
            for (i = 1; i <= NUMSAMPLES; i++) {
               LINE = LINE sprintf("\t%6i", AB[SAMPLE[ORDER[i]]][sprintf("%.1f", ab)] + 0)
            }
            print LINE
         }
      }' $VCFFILE > AB_perSample.hist
fi

# One more thing I can do is to estimate Hardy-Weinberg statistics for every
# population. The excess of homozygous individuals should deviate the populations
# further from HWE.

for pop in europaeus roumanicus concolor; do
   if [ ! -e $pop.hwe ] && [ ! -e $pop'_hwe.summary' ]; then
      vcftools --vcf $VCFFILE \
               --out $pop \
               --keep $pop.txt \
               --hardy &
   fi
done

wait

for  pop in europaeus roumanicus concolor; do
   if [ ! -e $pop'_hwe.summary' ]; then
      cut -f 3- $pop.hwe | grep -vP "nan|HWE" | gawk '{F[$0]++}END{for (f in F) print f "\t" F[f]}' | LC_ALL=C sort -grk 7 -k 3 > $pop'_hwe.summary'
   fi
done

if [ ! -e summary_hwe.txt ]; then
   echo "# These are numbers of sites polymorphic within a population that show significatn deficit" > summary_hwe.txt
   echo "# or excess of heterozygous individuals, or not, at 0.05 significance level, according to" >> summary_hwe.txt
   echo "# the tests performed by vcftools. Sites may be counted more than once, when polymorphic in" >> summary_hwe.txt
   echo "# more than one population." >> summary_hwe.txt
   gawk '($5 < 0.05){HET_DEF += $7}
      ($6 < 0.05){HET_EXC += $7}
      (($5 >= 0.05) && ($6 >= 0.05)){HWE += $7}
      END{
         TOTAL = HET_DEF + HET_EXC + HWE
         printf("Heterozygous deficit:%7i (%3.4f\%)\nHeterozyous excess:%9i (%3.4f\%)\nIn HW equilibrium:%10i (%3.4f\%)\n", HET_DEF, 100 * HET_DEF / TOTAL, HET_EXC, 100 * HET_EXC / TOTAL, HWE, 100 * HWE / TOTAL)
      }' *_hwe.summary >> summary_hwe.txt
fi

if [ ! -e summary_duplicate_reads.txt ]; then
   echo -e "#Library\tTotal_reads\tDuplicates\tPercentage" > summary_duplicate_reads.txt
   gawk -F "\t" '(($1 > 2) && ($1 <= 200)){
      DUPS[FILENAME] += $2
   }(NF == 3){
      TOTAL[FILENAME] = $3
   }END{
      for (f in DUPS) {
         LIBRARY = f
         gsub(/_dups.hist/, "", LIBRARY)
         print LIBRARY "\t" TOTAL[f] "\t" DUPS[f] "\t" 100 * DUPS[f] / TOTAL[f]
      }
   }' *.hist | sort -k 1,1 >> summary_duplicate_reads.txt
fi

# CONCLUSIONS
#
# From the bam files, I counted as duplicates the second reads that are aligned to the
# same position. I built the histograms of number of duplicates per positions, since it
# is very different to have 1000 duplicates from the same or from 500 different positions.
# I only count as problematic duplicates those that repeat a site more than twice and less
# than 201. The rationale is that adding one read to a locus will not bias much the allele
# frequency. And adding more than 200 duplicated reads to a locus will increase the coverage
# of all variant genotypes there above the threshold used to filter genotypes by coverage.
# The file "summary_duplicate_reads.txt" shows an important batch effect. The number of
# problematic duplicates drops with subsequent library preparations. The first 25 samples
# from 2016 have the highest levels of duplicates (between 5.8 and 7.8 %). The subsequent
# 25 samples from 2016 have 2.1 and 3.2 % duplicates. And all samples from 2018 have only
# between 0.4 and 0.7%.
#
# The allele balance (AB) in heterozygous sites is the proportion of each allele's observations.
# It should be around 50% of each allele. PCR duplicates increase the proportion of observations
# of one allele, biasing the AB at each site. The vcf file reports a measure of AB, defined as
# follows: "Allele balance at heterozygous sites: a number between 0 and 1 representing the
# ratio of reads showing the reference allele to all reads, considering only reads from individuals
# called as heterozygous". The average AB reported in the vcf file is 0.48 (see file 'AB.hist'),
# giving the false impression that PCR duplicates do not affect it. However, the allele that
# gets amplified by PCR duplicates is not always the same. The individual biases balance out
# when combining coverage data from several heterozygous individuals. The way to use AB as an
# indicator of PCR duplicates is per individual (see file 'AB_perSample.hist'. Then we see that
# the average proportion of observations of the least observed allele in heterozygous sites
# ranges from 0.21 to 0.36 among individuals. That's much less than the 50% expected without
# PCR duplicates.
#
# This suggests that several heterozygous sites are erroneously considered homozygous, because
# of undersampling of one allele. Granted that the allele missed is not always the same, so
# that the allele frequencies are still estimated accurately. But the genotype proportions and
# observed heterozygosities are biased.
#
# The Hardy-Weinberg statistics reported by vcftools (see 'summary_hwe.txt') show that 16% of
# sites polymorphic within a population show a significant (p value < 0.05) deficit of heterozygous
# genotypes. Only 0.7 % shows a significant excess.
