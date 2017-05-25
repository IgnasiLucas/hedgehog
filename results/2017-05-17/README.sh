#!/bin/bash
#
#				2017-05-17
#				----------
#
# Here I want to answer the question "Why there is not enough data to run
# most analyses?"

SABRELOG=/data/kristyna/hedgehog/data/runII_fastq/sabre.log
RUN_ONE=/data/kristyna/hedgehog/data/runI_fastq
BAMFILES=/data/kristyna/hedgehog/data/sorted_bam

# First, how much data is there? For the moment, we have the already demultiplexed
# reads from run I, and also the raw data from run II. Counting only first reads,
# From the log file of the demultiplexing of run II with sabre, I see the following:

SABRELOG=/data/kristyna/hedgehog/data/runII_fastq/sabre.log

gawk '(/FastQ records for barcode/){
   WITH += $6 / 2
}(/^FastQ records with no barcode/){
   WITHOUT += $7 / 2
}END{
   TOTAL = WITH + WITHOUT
   print "Reads_with_barcode\tReads_without_barcode\tTotal"
   printf "%u (%.2f %%)\t%u (%.2f %%)\t%u\n", WITH, 100 * WITH / TOTAL, WITHOUT, 100 * WITHOUT / TOTAL, TOTAL
}' $SABRELOG

# Output:
#
# Reads_with_barcode	Reads_without_barcode	Total
# 175901016 (67.47 %)	84827578 (32.53 %)	260728594
#
#
# From the demultiplexed fastq files of run I, I can count the lines and divide by
# 4 to count the number of first reads:

if [ ! -e read_counts.txt ]; then
   echo -e "Run\tSample\tOriginal_reads" > read_counts.txt
fi
for i in `ls -1 $RUN_ONE/*_1.fastq.gz`; do
   if ! grep -q `basename ${i:40} .fastq.gz` read_counts.txt; then
      zless $i | wc -l | gawk -v S=`basename ${i:40} .fastq.gz` '{print "run_I\t" S "\t" $1 / 4}' >> read_counts.txt
   fi
done

# The total number of reads in run I is 144033031. Among the two runs, we have 319934047
# reads.
#
# The distribution of reads per sample in the second run are given in the sabre log
# file.
if ! grep -q Er51_436 read_counts.txt; then
   gawk 'BEGIN{
      CODE["CGCGC"] = "Er51_436"
      CODE["CGGCG"] = "Er52_451"
      CODE["CGTAT"] = "Er53_ASR7"
      CODE["CTAGG"] = "Er54_AU1"
      CODE["CTCTT"] = "Er55_AU7"
      CODE["CTTCC"] = "Er56_AZ5"
      CODE["GAAGC"] = "Er57_COR4"
      CODE["GACTA"] = "Er58_FI7"
      CODE["GAGAT"] = "Er59_FR1"
      CODE["GATCG"] = "Er60_GR5"
      CODE["GCATT"] = "Er61_GR87"
      CODE["GCCGG"] = "Er62_GR95"
      CODE["GCGCC"] = "Er63_IR6"
      CODE["GCTAA"] = "Er64_IS1"
      CODE["GGAAG"] = "Er65_IS25"
      CODE["GGCCT"] = "Er66_IT3"
      CODE["GGGGA"] = "Er67_IT5"
      CODE["GGTTC"] = "Er68_PRT1B"
      CODE["GTACA"] = "Er69_R2"
      CODE["GTCAC"] = "Er70_RMN42"
      CODE["GTGTG"] = "Er71_SAR2"
      CODE["GTTGT"] = "Er72_SIE1B"
      CODE["TAATG"] = "Er73_SNG1"
      CODE["TACGT"] = "Er74_SP16"
      CODE["TAGCA"] = "Er75_TRC2A"
   }(/^FastQ records for barcode/){
      gsub(/:/,"",$5)
      print "run_II\t" CODE[$5] "\t" $6 / 2
   }' $SABRELOG >> read_counts.txt
fi

# To count the number of mapped reads, I use the sorted bam files
SAMPLE=(Er26_JUG4 Er27_SK32 Er28_112   Er29_122  Er30_79    Er31_453  Er32_183   Er33_211  Er34_197  Er35_209
        Er36_SK24 Er37_SK27 Er38_LB1   Er39_PL1  Er40_M2    Er41_GR36 Er42_GR35  Er43_SL7  Er44_VOJ1 Er45_BLG3
        Er46_RMN7 Er47_CR4  Er48_BH16  Er49_GR38 Er50_R3    Er51_436  Er52_451   Er53_ASR7 Er54_AU1  Er55_AU7
        Er56_AZ5  Er57_COR4 Er58_FI7   Er59_FR1  Er60_GR5   Er61_GR87 Er62_GR95  Er63_IR6  Er64_IS1  Er65_IS25
        Er66_IT3  Er67_IT5  Er68_PRT1B Er69_R2   Er70_RMN42 Er71_SAR2 Er72_SIE1B Er73_SNG1 Er74_SP16 Er75_TRC2A)

if [ ! -e mapped_counts.txt ]; then
   echo -e "Sample\tMapped" > mapped_counts.txt
   for i in ${SAMPLE[@]}; do
      samtools view -f 64 -F 260 $BAMFILES/$i.sortnc.bam | wc -l | gawk -v S=$i '{print S "\t" $1}' >> mapped_counts.txt
   done
fi
