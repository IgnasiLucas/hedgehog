#!/bin/bash
#
#				2019-02-27
#				----------
#
# I want to use snpEff to annotate the variants, in order to be able to run a
# McDonald-Kreitman test (or variation thereof) along the genome of E. roumanicus.
# To run this, you need python 2. I have prepared a conda environment specific
# for this folder with python 2.7 in spec-file.txt.

GZVCF=/data/kristyna/hedgehog/results_2018/23-02-2018/merged.vcf.gz

# SET UP
# ------
#
# The mapping of reads was originally done to the most recent version of the hedgehog
# reference genome, which is EriEur2.0. Thus, the names of contigs and coordinates of all
# variants are relative to that reference. An anotation of the genome is also available
# in the same coordinates at ftp://ftp.ncbi.nlm.nih.gov/genomes/Erinaceus_europaeus/GFF/ref_EriEur2.0_top_level.gff3.gz.
# But the program snpEff would import the annotation and genome sequence from a previous
# version, available in Ensembl. Thus, I need to build the EriEur2.0 database in snpEff.
# Before that, I need to add the "start_codon" and "stop_codon" features to the EriEur2.0
# gff3 file, because otherwise snpEff does not predict synonymous and non-synonymous variants.
# I found only one program that can add these features, which is called GAG. Below, I
# download and locally install both gag and snpEff to make the whole thing reproducible.
#
# This loop is a way to check for the existence of a file (a directory in this case),
# when I only know part of the name of that file (or directory).
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

if [ ! -d ../../bin/snpEff ]; then
   wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
   unzip snpEff_latest_core.zip -d ../../bin/
   rm snpEff_latest_core.zip
fi

if [ ! -d ../../bin/snpEff/data ]; then mkdir ../../bin/snpEff/data; fi
if [ ! -d ../../bin/snpEff/data/EriEur2.0 ]; then mkdir ../../bin/snpEff/data/EriEur2.0; fi
if [ ! -d ../../bin/snpEff/data/genomes ]; then mkdir ../../bin/snpEff/data/genomes; fi

# The names of scaffolds in the fasta file must coincide with those in the gff
if [ ! -e ../../bin/snpEff/data/genomes/EriEur2.0.fa ]; then
   wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Erinaceus_europaeus/CHR_Un/eeu_ref_EriEur2.0_chrUn.fa.gz
   gunzip -c eeu_ref_EriEur2.0_chrUn.fa.gz | gawk -v FS="|" '(/^>/){
      print ">" $4
   }(/^[^>]/){
      print $0
   }' > ../../bin/snpEff/data/genomes/EriEur2.0.fa
   rm eeu_ref_EriEur2.0_chrUn.fa.gz
fi

if [ ! -e ../../bin/snpEff/data/EriEur2.0/genes.gff ]; then
   if [ ! -e EriEur2.0.gff3 ]; then
      if [ ! -e ref_EriEur2.0_top_level.gff3.gz ]; then
         wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Erinaceus_europaeus/GFF/ref_EriEur2.0_top_level.gff3.gz
      fi
      gunzip -c ref_EriEur2.0_top_level.gff3.gz > EriEur2.0.gff3
      rm ref_EriEur2.0_top_level.gff3.gz
   fi
   if python --version 2>&1 | grep -q "2."; then
      python $GAG_DIR/gag.py --fasta ../../bin/snpEff/data/genomes/EriEur2.0.fa \
                             --gff EriEur2.0.gff3 \
                             --fix_start_stop \
                             --out start_stop_fixed

      mv start_stop_fixed/genome.gff ../../bin/snpEff/data/EriEur2.0/genes.gff
      mv start_stop_fixed/genome.stats ./
      rm -r start_stop_fixed
      rm EriEur2.0.gff3
   else
      echo "To run GAG, you need python 2."
      exit
   fi
fi

if ! grep -q EriEur2.0 ../../bin/snpEff/snpEff.config; then
   echo "EriEur2.0.genome : Erinaceus_europaeus" >> ../../bin/snpEff/snpEff.config
   echo "EriEur2.0.reference : ftp://ftp.ncbi.nlm.nih.gov/genomes/Erinaceus_europaeus/GFF/ref_EriEur2.0_top_level.gff3.gz" >> ../../bin/snpEff/snpEff.config
fi

# ANNOTATION
# ----------
if [ ! -e merged.annotated.vcf ]; then
   java -jar ../../bin/snpEff/snpEff.jar build -c ../../bin/snpEff/snpEff.config -gff3 -v EriEur2.0
   java -jar ../../bin/snpEff/snpEff.jar -c ../../bin/snpEff/snpEff.config -v EriEur2.0 $GZVCF > merged.annotated.vcf
fi

# SUMMARY
# -------

