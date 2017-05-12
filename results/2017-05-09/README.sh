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
   if [ ! -e ${SAMPLE[$i]}.fastq ] && [ ! -e whole.fastq ]; then
      vsearch --fastx_subsample $FASTQ_I_DIR/${SAMPLE[$i]}_1.fastq.gz --fastqout ${SAMPLE[$i]}.fastq --sample_size 1000 &
   fi
done

for i in `seq 25 49`; do
   if [ ! -e ${SAMPLE[$i]}.fastq ] && [ ! -e whole.fastq ]; then
      vsearch --fastx_subsample $FASTQ_II_DIR/${SAMPLE[$i]}_1.fastq.gz --fastqout ${SAMPLE[$i]}.fastq --sample_size 1000 &
   fi
done
wait

if [ ! -e whole.fastq ]; then
   cat Er*.fastq > whole.fastq
   rm Er*.fastq
fi

# I trim the reads with gawk because I notice some are already trimmed, and therefore
# I cannot just use the seqcrumb trim_edges.
if [ ! -e trimmed.fastq ]; then
   gawk '{$1 = substr($1,1,80); print}' whole.fastq > trimmed.fastq
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
   for j in whole trimmed; do
      for k in sensitive conservative; do
         if [ ! -d ${FOLDER[$i]}/$j.$k ]; then
            mkdir ${FOLDER[$i]}/$j.$k
         fi
         if [ ! -e ${FOLDER[$i]}/$j.$k/output.log ]; then
            ./align.sh ${FOLDER[$i]} ${EXE[$i]} $j.fastq $k ${FOLDER[$i]}/$j.$k 2> ${FOLDER[$i]}/$j.$k/output.log
#         else
#            echo "${FOLDER[$i]}/$j.$k/output.log already present."
#            exit
         fi
      done
   done
done

if [ ! -e summary.txt ]; then
   echo -e "Aligner\tInput\tSetting\tUnique\tMultiple\tUnmapped" > summary.txt
   for i in whole trimmed; do
      for j in conservative sensitive; do
         gawk -v INPUT=$i -v SETTING=$j '(/reads processed:/){
            TOTAL = $4
         }(/reads with at least one reported alignment:/){
            UNIQUE = $9
         }(/reads that failed to align:/){
            UNMAPPED = $7
         }(/reads with alignments suppressed due to -m:/){
            MULTIPLE = $9
         }END{
            printf "Bowtie1\t%s\t%s\t%.2f%%\t%.2f%%\t%.2f%%\n", INPUT, SETTING, 100 * (UNIQUE + 0) / TOTAL, 100 * (MULTIPLE + 0) / TOTAL, 100 * (UNMAPPED + 0) / TOTAL
         }' bowtie1/$i.$j/output.log >> summary.txt
      done
   done
   for i in whole trimmed; do
      for j in conservative sensitive; do
         gawk -v INPUT=$i -v SETTING=$j '(/reads; of these:/){
            TOTAL = $1
         }(/aligned exactly 1 time/){
            UNIQUE = $1
         }(/aligned 0 times/){
            UNMAPPED = $1
         }(/aligned >1 times/){
            MULTIPLE = $1
         }END{
            printf "Bowtie2\t%s\t%s\t%.2f%%\t%.2f%%\t%.2f%%\n", INPUT, SETTING, 100 * (UNIQUE + 0) / TOTAL, 100 * (MULTIPLE + 0) / TOTAL, 100 * (UNMAPPED + 0) / TOTAL
         }' bowtie2/$i.$j/output.log >> summary.txt
      done
   done
fi

# CONCLUSIONS
# ===========
#
# Summary table:
#
# --------------------------------------------------------------
# Aligner	Input	Setting 	Unique	Multip.	Unmapped
# --------------------------------------------------------------
# Bowtie1	whole	conservative	46.41%	6.57%	47.02%
# Bowtie1	whole	sensitive	59.73%	0.00%*	40.27%
# Bowtie1	trimmed	conservative	48.96%	13.22%	37.81%
# Bowtie1	trimmed	sensitive	69.10%	0.00%*	30.90%
# Bowtie2	whole	conservative	47.22%	34.32%	18.46%
# Bowtie2	whole	sensitive	45.51%	32.63%	21.86%
# Bowtie2	trimmed	conservative	36.70%	47.03%	16.26%
# Bowtie2	trimmed	sensitive	36.36%	45.98%	17.66%
# --------------------------------------------------------------
#
# * In the sensitive setting of bowtie1 the -m options, which removes ambiguously mapped reads
# was not applied. Thus, there are no multiply mapped reads identified, but they are actually
# counted together with the uniquely mapped.
#
# In general, the lower percentage of unmapped reads with Bowtie2 confirms that Bowtie2 is
# more sensitive than Bowtie1. However, it is clear that the higher sensitivity only helps
# finding secondary mappings to many reads.
#
# The current mappings, made with Bowtie1 produced an amount of uniquely mapped reads similar
# to the best results tested here. Thus, it is not necessary to try to improve the mapping.
# Instead of increasing coverage, we would increase the noise caused by ambiguously mapped
# reads.
