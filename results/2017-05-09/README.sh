#!/bin/bash
#
#				2017-05-09
#				----------
#
# I'm going to test some alignment strategies of the first reads, in order to decide
# if it is necessary to align again the whole dataset. I depart from the fastq files
# from run I and run II, which are already demultiplexed in the data folder.
#
# -----------------------------------------------------
# SET PATHS:
# -----------------------------------------------------

  FASTQ_I_DIR=/data/kristyna/hedgehog/data/runI_fastq
 FASTQ_II_DIR=/data/kristyna/hedgehog/data/runII_fastq
    REFERENCE=/data/kristyna/hedgehog/data/reference/GCF_000296755.1_EriEur2.0_genomic.fna.gz
       BOWTIE=/home/joiglu/miniconda2/bin/bowtie
 BOWTIE_BUILD=/home/joiglu/miniconda2/bin/bowtie-build
      BOWTIE2=/home/joiglu/miniconda2/bin/bowtie2
BOWTIE2_BUILD=/home/joiglu/miniconda2/bin/bowtie2-build

# -----------------------------------------------------


SAMPLE=(Er26_JUG4 Er27_SK32 Er28_112   Er29_122  Er30_79    Er31_453  Er32_183   Er33_211  Er34_197  Er35_209
        Er36_SK24 Er37_SK27 Er38_LB1   Er39_PL1  Er40_M2    Er41_GR36 Er42_GR35  Er43_SL7  Er44_VOJ1 Er45_BLG3
        Er46_RMN7 Er47_CR4  Er48_BH16  Er49_GR38 Er50_R3    Er51_436  Er52_451   Er53_ASR7 Er54_AU1  Er55_AU7
        Er56_AZ5  Er57_COR4 Er58_FI7   Er59_FR1  Er60_GR5   Er61_GR87 Er62_GR95  Er63_IR6  Er64_IS1  Er65_IS25
        Er66_IT3  Er67_IT5  Er68_PRT1B Er69_R2   Er70_RMN42 Er71_SAR2 Er72_SIE1B Er73_SNG1 Er74_SP16 Er75_TRC2A)


for i in `seq 0 24`; do
   if [ ! -e ${SAMPLE[$i]}.fastq ] && [ ! -e read1.fastq ]; then
      vsearch --fastx_subsample $FASTQ_I_DIR/${SAMPLE[$i]}_1.fastq.gz --fastqout ${SAMPLE[$i]}.fastq --sample_size 1000 &
   fi
done

for i in `seq 25 49`; do
   if [ ! -e ${SAMPLE[$i]}.fastq ] && [ ! -e read1.fastq ]; then
      vsearch --fastx_subsample $FASTQ_II_DIR/runII_fastq/${SAMPLE[$i]}_1.fastq.gz --fastqout ${SAMPLE[$i]}.fastq --sample_size 1000 &
   fi
done
wait

if [ ! -e read1.fastq ]; then
   cat Er*.fastq > read1.fastq
   rm Er*.fastq
fi

# I trim the reads with gawk because I notice some are already trimmed, and therefore
# I cannot just use the seqcrumb trim_edges.
if [ ! -e trimmed.fastq ]; then
   gawk '{$1 = substr($1,1,80); print}' read1.fastq > trimmed.fastq
fi

if [ ! -e bowtie1 ]; then mkdir bowtie1; fi
if [ ! -e bowtie2 ]; then mkdir bowtie2; fi


if [ ! -e bowtie1/erinaceus.1.ebwt ]; then
   # Apparently, neither bowtie nor bowtie2 index compressed fasta files.
   if [ ! -e z1.fa ]; then
      if [ ! -e z1.fa.gz ]; then
         cp $REFERENCE ./z1.fa.gz
      fi
      gunzip z1.fa.gz
   fi
   $BOWTIE_BUILD z1.fa bowtie1/erinaceus
fi

if [ ! -e bowtie2/erinaceus.1.bt2 ]; then
   if [ ! -e z1.fa ]; then
      if [ ! -e z1.fa.gz ]; then
         cp $REFERENCE ./z1.fa.gz
      fi
      gunzip z1.fa.gz
   fi
   $BOWTIE2_BUILD z1.fa bowtie2/erinaceus
   rm z1.fa
fi

EXE=($BOWTIE $BOWTIE2)
FOLDER=(bowtie1 bowtie2)
for i in 0 1; do
   for j in read1 trimmed; do
      for k in sensitive conservative; do
         if [ ! -d ${FOLDER[$i]}/$j.$k ]; then
            mkdir ${FOLDER[$i]}/$j.$k
         fi
         align.sh ${FOLDER[$i]} ${EXE[$i]} $j.fastq $k ${FOLDER[$i]}/$j.$k
      done
   done
done
