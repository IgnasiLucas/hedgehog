#!/bin/bash
#
# I need to clarify how many reads were produced, aligned, and nuclear in each
# sequencing run. The first two sequencing runs, both from 2016, targeted different
# individuals, 50 in all, and will be considered together. These are the paths to
# the data.
#
# FASTQ2016 include runs I and II, and they have already been demultiplexed.
FASTQ2016=/data/kristyna/hedgehog/data_2016/fastq

# BAM2016 were mapped with bowtie2. There is a set of alignments with bowtie as
# well. But I think we used the bowtie2 alignments. However, I do not see the
# nuclear subset of those alignments. The log files from bowtie2 are available,
# and give the number of reads mapped.
BAM2016=/data/kristyna/hedgehog/data_2016/bam

# FASTQ2018 includes demultiplexed reads from the third run. It includes the
# log file from sabre, which already counts how many reads there should be in the
# individual files. The raw data is also available.
FASTQ2018=/data/kristyna/hedgehog/data_2018/fastq

# BAM2018 includes the bam files from bowtie2, and also the log files.
BAM2018=/data/kristyna/hedgehog/data_2018/bam

# MERGED includes merged bam files from both sequencing seasons, and also the
# selection of nuclear reads (*_msnc.bam).
MERGED=/data/kristyna/hedgehog/data_2018/merged_bam

SAMPLE=(Er26_JUG4 Er27_SK32 Er28_112   Er29_122  Er30_79    Er31_453  Er32_183   Er33_211  Er34_197  Er35_209
        Er36_SK24 Er37_SK27 Er38_LB1   Er39_PL1  Er40_M2    Er41_GR36 Er42_GR35  Er43_SL7  Er44_VOJ1 Er45_BLG3
        Er46_RMN7 Er47_CR4  Er48_BH16  Er49_GR38 Er50_R3    Er51_436  Er52_451   Er53_ASR7 Er54_AU1  Er55_AU7
        Er56_AZ5  Er57_COR4 Er58_FI7   Er59_FR1  Er60_GR5   Er61_GR87 Er62_GR95  Er63_IR6  Er64_IS1  Er65_IS25
        Er66_IT3  Er67_IT5  Er68_PRT1B Er69_R2   Er70_RMN42 Er71_SAR2 Er72_SIE1B Er73_SNG1 Er74_SP16 Er75_TRC2A)


# Count the number of reads in the fastq file:
READ1_2016=$( gunzip -c $FASTQ2016/$1"_1.fastq.gz" | wc -l | gawk '{print $1 / 4}' )
READ2_2016=$( gunzip -c $FASTQ2016/$1"_2.fastq.gz" | wc -l | gawk '{print $1 / 4}' )
READ1_2018=$( gunzip -c $FASTQ2018/$1"_1.fastq.gz" | wc -l | gawk '{print $1 / 4}' )
READ2_2018=$( gunzip -c $FASTQ2018/$1"_2.fastq.gz" | wc -l | gawk '{print $1 / 4}' )

# Count the number of read pairs reported in bowtie2 log files as total number of reads:
LOGREADS_2016=$( gawk '(/reads; of these:/){print $1}' $BAM2016/$1"_bowtie2_2016.log" )
LOGREADS_2018=$( gawk '(/reads; of these:/){print $1}' $BAM2018/$1"_bowtie2_2018.log" )

# Count the number of reads successfully mapped, either concordantly or not. I am using
# two different ways to count them. From the bam file, I count primary mappings of first
# reads. From the log file, I count reads (pairs) mapped concordantly once or more times,
# and discordantly.
MAPPED1_2016=$( samtools view -f 64 -F 260 $BAM2016/$1"_2016.bam" | wc -l )
MAPPED1_2018=$( samtools view -f 64 -F 260 $BAM2018/$1"_2018.bam" | wc -l )
LOGMAPPED_2016=$( gawk '(/1 time/){S += $1}END{print S}' $BAM2016/$1"_bowtie2_2016.log" )
LOGMAPPED_2018=$( gawk '(/1 time/){S += $1}END{print S}' $BAM2018/$1"_bowtie2_2018.log" )

# Count the merged mapped reads, and how many of them are nuclear.
MAPPED_BOTH=$( samtools view -f 64 -F 260 $MERGED/$1"_merged.bam" | wc -l )
MAPPED_SORT=$( samtools view -f 64 -F 260 $MERGED/$1"_merged_sorted.bam" | wc -l )
 MAPPED_NUC=$( samtools view -f 64 -F 260 $MERGED/$1"_msnc.bam" | wc -l )

echo -e "$1\t$READ1_2016\t$READ2_2016\t$LOGREADS_2016\t$READ1_2018\t$READ2_2018\t$LOGREADS_2018\t$MAPPED1_2016\t$LOGMAPPED_2016\t$MAPPED1_2018\t$LOGMAPPED_2018\t$MAPPED_BOTH\t$MAPPED_SORT\t$MAPPED_NUC"
