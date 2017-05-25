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
   fi
done

if [ ! -e ref.fa ]; then
   gunzip -c $REFERENCE > ref.fa
fi

if [ ! -e hhog.vcf ]; then
   freebayes -f ref.fa \
      --vcf hhog.vcf \
      --populations populations.txt \
      --theta 0.01 \
      --standard-filters \
      --max-coverage 900 \
      --genotype-qualities \
      ./*.bam
fi
