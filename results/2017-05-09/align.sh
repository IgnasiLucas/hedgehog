#!/bin/bash
#
# This is called by README.sh like this:
#
#    align.sh ${FOLDER[$i]} ${EXE[$i]} $j.fastq $k ${FOLDER[$i]}/$j.$k
#
# Where arguments are:
#
#          ${FOLDER[$i]}:   either bowtie1 or bowtie2, where the indices are and results are sent.
#             ${EXE[$i]}:   the executable, either bowtie or bowtie2, with whole path.
#               $j.fastq:   the input fastq file, either trimmed or not.
#                     $k:   either "sensitive" or "conservative".
#    ${FOLDER[$i]}/$j.$k:   the subfolder where to sent the output.
#

if [ $1 == "-h" ]; then
   echo "./align.sh <index folder> <path to bowtie or bowtie2> <input fastq> {sensitive | conservative} <output dir>"
   exit
fi

if [ $1 == bowtie1 ]; then
   if [ $4 == sensitive ]; then
      if [ ! -e bowtie1/trimmed.sensitive/output.bam ]; then
         # I request the best alignment among those with no more than 3 mismatches in the first 60
         # bases, and with a total base quality sum of mismatching bases not higher than 80.
         $2 -n 3 -l 60 -e 80 --best -p 20 --nomaqround -t -S bowtie1/erinaceus $3 | samtools view -Sb - > $5/output.bam
      fi
      if [ ! -e bowtie1/whole.sensitive/output.bam ]; then
         $2 -n 3 -l 60 -e 80 --best -p 20 --nomaqround -t -S bowtie1/erinaceus $3 | samtools view -Sb - > $5/output.bam
      fi
   else
      if [ $4 == conservative ]; then
         if [ ! -e bowtie1/trimmed.conservative/output.bam ]; then
            # In the conservative setting, I request no more than 1 mismatch in the first 50 bases,
            # with a total base quality sum of mismatching bases not higher than 60. I also dismiss
            # ambiguous mappings.
            $2 -n 1 -l 50 -e 60 -m 1 -p 20 --nomaqround -t -S bowtie1/erinaceus $3 | samtools view -Sb - > $5/output.bam
         fi
         if [ ! -e bowtie1/whole.conservative/output.bam ]; then
            $2 -n 1 -l 50 -e 60 -m 1 -p 20 --nomaqround -t -S bowtie1/erinaceus $3 | samtools view -Sb - > $5/output.bam
         fi
      fi
   fi
else
   if [ $1 == bowtie2 ]; then
      if [ $4 == sensitive ]; then
         if [ ! -e bowtie2/trimmed.sensitive/output.bam ]; then
            # I request up to 1 mismatch in 20 bp long seeds. This is between the sensitive and
            # the very sensitive preset options.
            $2 -N 1 -L 20 -D 25 -R 2 -i S,1,1 -t --no-unal -p 20 -x bowtie2/erinaceus -U $3 | samtools view -Sb - > $5/output.bam
         fi
         if [ ! -e bowtie2/whole.sensitive/output.bam ]; then
            $2 -N 1 -L 20 -D 25 -R 2 -i S,1,1 -t --no-unal -p 20 -x bowtie2/erinaceus -U $3 | samtools view -Sb - > $5/output.bam
         fi
      else
         if [ $4 == conservative ]; then
            if [ ! -e bowtie2/trimmed.conservative/output.bam ]; then
               # Zero mismatches in 22-bp long seeds.
               $2 -N 0 -L 22 -D 10 -R 2 -i S,1,1.15 -t --no-unal -p 20 -x bowtie2/erinaceus -U $3 | samtools view -Sb - > $5/output.bam
            fi
            if [ ! -e bowtie2/whole.conservative/output.bam ]; then
               $2 -N 0 -L 22 -D 10 -R 2 -i S,1,1.15 -t --no-unal -p 20 -x bowtie2/erinaceus -U $3 | samtools view -Sb - > $5/output.bam
            fi
         fi
      fi
   fi
fi
