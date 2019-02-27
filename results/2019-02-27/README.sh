#!/bin/bash
#
#				2019-02-27
#				----------
#
# I want to use snpEff to annotate the variants, in order to be able to run a
# McDonald-Kreitman test (or variation thereof) along the genome of E. roumanicus.
#
# The databases preloaded in snpEff already included a hedgehog genome annotation.
# Unfortunately, it was the previous version of the reference genome that was used
# in snpEff, and I had to build a new database in snpEff, following the instructions
# in http://snpeff.sourceforge.net/SnpEff_manual.html#databases. I apologize for
# not keeping records on how I did that.
#
GZVCF=/data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz

# This loop is a way to check for the existence of a file (a directory in this case),
# when I only know part of the name of that file (or directory). I need to GAG program
# to add codon_start and codon_stop features to the gff.
for f in ../../bin/genomeannotation-GAG-*; do
   if [ -d $f ]; then
      GAG_DIR=$f
   else
      wget https://github.com/genomeannotation/GAG/zipball/master
      unzip master -d ../../bin/
      rm master
      GAG_DIR=$(find ../../bin -maxdepth 1 -type d -name 'genomeannotation*' | head -n 1)
   fi
   break
done

if [ ! -e genes.gff ]; then
   if [ ! -e EriEur2.0.gff3 ]; then
      if [ ! -e ref_EriEur2.0_top_level.gff3.gz ]; then
         wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Erinaceus_europaeus/GFF/ref_EriEur2.0_top_level.gff3.gz
      fi
      gunzip -c ref_EriEur2.0_top_level.gff3.gz > EriEur2.0.gff3
      rm ref_EriEur2.0_top_level.gff3.gz
   fi
   # The names of scaffolds in the fasta file must coincide with those in the gff
   if [ ! -e reference.fa ]; then
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Erinaceus_europaeus/CHR_Un/eeu_ref_EriEur2.0_chrUn.fa.gz
      gunzip -c eeu_ref_EriEur2.0_chrUn.fa.gz | gawk -v FS="|" '(/^>/){
         print ">" $4
      }(/^[^>]/){
         print $0
      }' > reference.fa
      rm eeu_ref_EriEur2.0_chrUn.fa.gz
   fi
   if python --version 2>&1 | grep -q "2."; then
      python $GAG_DIR/gag.py --fasta reference.fa --gff EriEur2.0.gff3 --fix_start_stop --out start_stop_fixed
      mv start_stop_fixed/genome.gff ./genes.gff
   else
      echo "To run GAG, you need python 2."
      exit
   fi
fi

if [ ! -d ../../bin/snpEff ]; then
   wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
   unzip snpEff_latest_core.zip -d ../../bin/
   rm snpEff_latest_core.zip
fi

if ! grep -q EriEur2.0 ../../bin/snpEff/snpEff.config; then
   echo "EriEur2.0.genome : Erinaceus_europaeus" >> ../../bin/snpEff/snpEff.config
   echo "EriEur2.0.reference : ftp://ftp.ncbi.nlm.nih.gov/genomes/Erinaceus_europaeus/GFF/ref_EriEur2.0_top_level.gff3.gz" >> ../../bin/snpEff/snpEff.config
fi

if [ ! -d ../../bin/snpEff/data ]; then mkdir ../../bin/snpEff/data; fi
if [ ! -d ../../bin/snpEff/data/EriEur2.0 ]; then mkdir ../../bin/snpEff/data/EriEur2.0; fi
if [ ! -e ../../bin/snpEff/data/EriEur2.0/genes.gff ]; then ln $(pwd)/genes.gff ../../bin/snpEff/data/EriEur2.0/genes.gff
if [ ! -d ../../bin/snpEff/data/genomes ]; then mkdir ../../bin/snpEff/data/genome; fi
if [ ! -e ../../bin/snpEff/data/genomes/EriEur2.0.fa ]; then ln $(pwd)/reference.fa ../../bin/snpEff/data/genomes/EriEur2.0.fa

if [ ! -e merged.annotated.vcf ]; then
   java -jar ../../bin/snpEff/snpEff.jar build -c ../../bin/snpEff/snpEff.config -gff3 -v EriEur2.0
   java -jar ../../bin/snpEff/snpEff.jar -c ../../bin/snpEff/snpEff.config -v EriEur2.0 $GZVCF > merged.annotated.vcf
fi
