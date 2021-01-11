#!/bin/bash
#
#				2020-06-15
#				==========
#
DATA="../2018-04-13/both_outfiles/both.alleles.loci"
POPMAP="../../data/populations.txt"

# First, I want to create fasta files for all loci. The biostrings package can
# deal with extended IUPAC symbols, representing heterozygous sites. However,
# the information is not easy to translate into correct number of differences
# among sequences. It is better to use the 'both.alleles.loci' file, with two
# sequences per individual and locus.
#
# There are 164546 alignments or loci, but they are numbered not
# consecutively, up to 834620. Using grep, I see in the summary line after every
# alignment that there are 140558 alignments with at least 1 and at most 10 variable
# sites. I will keep only those.

if [ ! -d fasta ]; then mkdir fasta; fi
if [ ! -e fasta/L834620.fa ]; then
   gawk '(/^[^\/]/){
      SEQUENCE[$1] = $2
   }(/^\//){
      if ((NF > 2) && (NF <= 12)) {
         gsub(/*|-|\|/, "", $NF)
         FASTAFILE = sprintf("fasta/L%06u.fa", $NF)
         for (ind in SEQUENCE) {
            print ">" ind "\n" SEQUENCE[ind] >>FASTAFILE
         }
      }
      delete SEQUENCE
   }' $DATA
fi

if [ ! -e mka.html ]; then
   R -q --save -e "rmarkdown::render('mka.Rmd', output_file='mka.html')"
fi
