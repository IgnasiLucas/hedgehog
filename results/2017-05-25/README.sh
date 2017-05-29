#!/bin/bash
#
#				2017-05-25
#				----------
#
# Kristýna used bowtie (version 1) to align trimmed reads to the reference
# genome. The sorted bam files of the nuclear alignments are available, but
# they miss the read group information. Here, I will use picard tools to
# add read group information and then I will used freebayes to call variants.

BAMDIR=/data/kristyna/hedgehog/data/sorted_bam
REFERENCE=/data/kristyna/hedgehog/reference/GCF_000296755.1_EriEur2.0_genomic.fna.gz
PICARD=/usr/local/bin/picard.jar
SAMPLE=(Er26_JUG4 Er27_SK32 Er28_112   Er29_122  Er30_79    Er31_453  Er32_183   Er33_211  Er34_197  Er35_209
        Er36_SK24 Er37_SK27 Er38_LB1   Er39_PL1  Er40_M2    Er41_GR36 Er42_GR35  Er43_SL7  Er44_VOJ1 Er45_BLG3
        Er46_RMN7 Er47_CR4  Er48_BH16  Er49_GR38 Er50_R3    Er51_436  Er52_451   Er53_ASR7 Er54_AU1  Er55_AU7
        Er56_AZ5  Er57_COR4 Er58_FI7   Er59_FR1  Er60_GR5   Er61_GR87 Er62_GR95  Er63_IR6  Er64_IS1  Er65_IS25
        Er66_IT3  Er67_IT5  Er68_PRT1B Er69_R2   Er70_RMN42 Er71_SAR2 Er72_SIE1B Er73_SNG1 Er74_SP16 Er75_TRC2A)

for i in ${SAMPLE[@]}; do
   if [ ! -e $i.bam ]; then
      java -jar $PICARD AddOrReplaceReadGroups \
         I=$BAMDIR/$i.sortnc.bam \
         O=z1.bam \
         RGID=$i \
         RGLB=$i \
         RGPL=illumina \
         RGPU=unit1 \
         RGSM=$i
      # I remove unmapped and second reads:
      samtools view -b -F 132 z1.bam > $i.bam
      rm z1.bam
   fi
done

if [ ! -e ref.fa ]; then
   gunzip -c $REFERENCE > ref.fa
fi

if [ ! -e hhog.vcf ]; then
   if [ ! -e z1.vcf ]; then
      freebayes -f ref.fa \
         --vcf z1.vcf \
         --populations populations.txt \
         --theta 0.01 \
         --standard-filters \
         --max-coverage 900 \
         --genotype-qualities \
         ./*.bam
   fi
   gawk -f add_flag.awk z1.vcf > hhog.vcf
   rm z1.vcf
fi

# Now that I have the vcf file, with binary presence flags, I should get some
# statistics that guide the filtering of variable sites. First, the histogram of the
# number of samples with data per locus.

if [ ! -e SamplesWithData.txt ]; then
   echo -e "#withData\tFrequency" > SamplesWithData.txt
   gawk '(/^[^#]/){
      split($8, INFO, /;/)
      for (i in INFO) {
         if (INFO[i] ~ /^NS=/) {
            F[substr(INFO[i],4)]++
         }
      }
   }END{
      for (f = 1; f <= 50; f++) print f "\t" F[f] + 0
   }' hhog.vcf | sort -nk 1,1 >> SamplesWithData.txt
fi

if [ ! -e SamplesWithData.png ]; then
   gnuplot < SamplesWithData.gnp
fi

# Below I use some commands to find the binary expression of the set of
# samples from each population. Then I will use that to filter the hhog.vcf
# according to how many samples from each population have data.

# This is the order in which the samples are in hhog.vcf:
head -n 10000 hhog.vcf | grep "#CHROM" | sed 's/\t/\n/g' | grep "^Er" | gawk '{print NR "\t" $1}' | sort -k 2,2 > z1

# Ordering populations.txt by the same criteria, will produce a parallele list of
# sample names, with the assigned population:

sort -k 1,1 populations.txt > z2

# now, we paste, order by the original positions in hhog.vcf, and keep the list of
# populations:
paste z1 z2 | sort -nk 1,1 | gawk '(/europaeus/){
      EUROPAEUS += 2^(NR - 1)
   }(/romanicus/){
      ROMANICUS += 2^(NR - 1)
   }(/concolor/){
      CONCOLOR  += 2^(NR - 1)
   }(/Atelerix/){
      ATELERIX  += 2^(NR - 1)
   }(/Hemiechinus/){
      HEMIECHINUS += 2^(NR - 1)
   }(/hybrid/){
      HYBRID += 2^(NR - 1)
   }END{
      print "Europaeus\t" EUROPAEUS
      print "Romanicus\t" ROMANICUS
      print "Concolor\t" CONCOLOR
      print "Atelerix\t" ATELERIX
      print "HemiEchin\t" HEMIECHINUS
      print "Hybrid\t" HYBRID
   }'

rm z1 z2

# I include this information in the script filtervcf.awk and I require at least 2
# europaeus, 2 romanicus, and 1 concolor in the output. I get a reduction of the
# vcf file size of about 50%.

if [ ! -e hhog_min2each.vcf ]; then
   gawk -f filtervcf.awk hhog.vcf > hhog_min2each.vcf
fi

