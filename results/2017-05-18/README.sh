#!/bin/bash
#
#				2017-05-18
#				----------
#
# Let's run ipyrad. I just updated it to version 0.6.19.

RUNONE=/data/kristyna/hedgehog/data/runI_fastq
RUNTWO=/data/kristyna/hedgehog/data/runII_fastq
REFERENCE=/data/kristyna/hedgehog/data/reference/GCF_000296755.1_EriEur2.0_genomic.fna.gz
# ipyrad needs write persion in the reference folder, to do the bwa index. In addition,
# the reference must be either uncompressed fasta or bgzip-compressed. Not gzip-compressed.

if [ ! -e reference.fasta ]; then
   cp $REFERENCE reference.fasta.gz
   gunzip reference.fasta.gz
fi

FASTQDIR=`pwd | sed "s/results/data/"`"/fastq"
if [ ! -e `pwd | sed "s/results/data/"` ]; then mkdir `pwd | sed "s/results/data/"`; fi
if [ ! -e $FASTQDIR ]; then mkdir $FASTQDIR; fi
SAMPLE=(Er26_JUG4 Er27_SK32 Er28_112   Er29_122  Er30_79    Er31_453  Er32_183   Er33_211  Er34_197  Er35_209
        Er36_SK24 Er37_SK27 Er38_LB1   Er39_PL1  Er40_M2    Er41_GR36 Er42_GR35  Er43_SL7  Er44_VOJ1 Er45_BLG3
        Er46_RMN7 Er47_CR4  Er48_BH16  Er49_GR38 Er50_R3    Er51_436  Er52_451   Er53_ASR7 Er54_AU1  Er55_AU7
        Er56_AZ5  Er57_COR4 Er58_FI7   Er59_FR1  Er60_GR5   Er61_GR87 Er62_GR95  Er63_IR6  Er64_IS1  Er65_IS25
        Er66_IT3  Er67_IT5  Er68_PRT1B Er69_R2   Er70_RMN42 Er71_SAR2 Er72_SIE1B Er73_SNG1 Er74_SP16 Er75_TRC2A)

for i in `seq 0 24`; do
   if [ ! -e $FASTQDIR/${SAMPLE[$i]}.fastq.gz ]; then
      ln -s $RUNONE/${SAMPLE[$i]}'_1.fastq.gz' $FASTQDIR/${SAMPLE[$i]}.fastq.gz
   fi
done
for i in `seq 25 49`; do
   if [ ! -e $FASTQDIR/${SAMPLE[$i]}.fastq.gz ]; then
      ln -s $RUNTWO/${SAMPLE[$i]}'_1.fastq.gz' $FASTQDIR/${SAMPLE[$i]}.fastq.gz
   fi
done

if [ ! -e popmap.txt ]; then
   echo "Populations file not found."
   exit
fi

if [ ! -e params-hhog.txt ]; then
   ipyrad -n hhog
   echo >> params-hhog.txt
#                      ------- ipyrad params file (v.0.6.19)-------------------------------------------
#  sed -i  "/## \[0\]/c hhog                           ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps"         params-hhog.txt
#  sed -i  "/## \[1\]/c ./                             ## [1] [project_dir]: Project dir (made in curdir if not present)"                                 params-hhog.txt
#  sed -i  "/## \[2\]/c                                ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files"                            params-hhog.txt
#  sed -i  "/## \[3\]/c                                ## [3] [barcodes_path]: Location of barcodes file"                                                 params-hhog.txt
   sed -i  "/## \[4\]/c $FASTQDIR/*.fastq.gz           ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files"                          params-hhog.txt
   sed -i  "/## \[5\]/c denovo+reference               ## [5] [assembly_method]: Assembly method (denovo, reference, denovo+reference, denovo-reference)" params-hhog.txt
   sed -i  "/## \[6\]/c ./reference.fasta              ## [6] [reference_sequence]: Location of reference sequence file"                                  params-hhog.txt
   sed -i  "/## \[7\]/c rad                            ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc."                                     params-hhog.txt
   sed -i  "/## \[8\]/c TGCAGG                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)"                       params-hhog.txt
   sed -i  "/## \[9\]/c 5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read"                          params-hhog.txt
   sed -i "/## \[10\]/c 33                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)"             params-hhog.txt
   sed -i "/## \[11\]/c 6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling"                            params-hhog.txt
   sed -i "/## \[12\]/c 5                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling"                              params-hhog.txt
   sed -i "/## \[13\]/c 150                            ## [13] [maxdepth]: Max cluster depth within samples"                                              params-hhog.txt
   sed -i "/## \[14\]/c 0.85                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly"                              params-hhog.txt
#  sed -i "/## \[15\]/c 0                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes"                    params-hhog.txt
   sed -i "/## \[16\]/c 0                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)"                          params-hhog.txt
#  sed -i "/## \[17\]/c 35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim"                             params-hhog.txt
   sed -i "/## \[18\]/c 2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences"                        params-hhog.txt
   sed -i "/## \[19\]/c 5, 5                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)"                          params-hhog.txt
   sed -i "/## \[20\]/c 8, 8                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)"                            params-hhog.txt
   sed -i "/## \[21\]/c 4                              ## [21] [min_samples_locus]: Min # samples per locus for output"                                   params-hhog.txt
   sed -i "/## \[22\]/c 20, 20                         ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)"                                           params-hhog.txt
   sed -i "/## \[23\]/c 3, 3                           ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)"                                    params-hhog.txt
   sed -i "/## \[24\]/c 0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)"                        params-hhog.txt
   sed -i "/## \[25\]/c 0, -5, 0, 0                    ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)"                         params-hhog.txt
   sed -i "/## \[26\]/c 0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)"                             params-hhog.txt
   sed -i "/## \[27\]/c *                              ## [27] [output_formats]: Output formats (see docs)"                                               params-hhog.txt
   sed -i "/## \[28\]/c popmap.txt                     ## [28] [pop_assign_file]: Path to population assignment file"                                     params-hhog.txt
fi

for STEP in 1 2 3 4 5 6 7; do
   if ipyrad -p params-hhog.txt -r | gawk -v STEP=$STEP '((/^step/) && ($2 == STEP ":")){print}' | tee z1 | grep -q None; then
      ipyrad -p params-hhog.txt -s $STEP -c 24 1> step$STEP.log 2> step$STEP.err
   else
      cat z1
      rm z1
   fi
done
