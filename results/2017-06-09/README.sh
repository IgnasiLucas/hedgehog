#!/bin/bash
#
#				2017-06-09
#				----------
#
# Loci with more than two parsimony informative sites are easy to filter
# from the alleles.loci output of ipyrad. I should output each alignment
# to a different file. I think it is safe to skip the second allele of an
# individual if it is identical to the first one. I do not want too many
# sequences to analyse, and most such cases must indicate undersampling,
# rather than real homozygosity.
#
# I cannot use the nexus output file from ipyrad, because it concatenates
# all loci, assuming that all loci must be sampled from the same set of
# individuals. Actually, it is interesting to note that when keeping loci independent
# it does not matter if the chromosomes sampled at each locus are different.
#
# It is probably easy to parse the output alignments directly from within the
# API of ipyrad, using python commands to filter and output the alignments in
# different nexus files. But the API is not well documented, and I have no
# idea how to do that.  It's easy enough to use awk.

LOCI=../2017-05-18/min5_outfiles/min5.alleles.loci
JMODELTEST_HOME=/usr/local/jmodeltest2/dist/
POPMAP=../2017-05-18/popmap.txt

if [ ! -e popmap.txt ]; then
   gawk 'BEGIN{
      SPECIES["Atelerix"]    = "Atelerix_sp"
      SPECIES["Hemiechinus"] = "Hemiechinus_sp"
      SPECIES["concolor"]    = "Erinaceus_concolor"
      SPECIES["romanicus"]   = "Erinaceus_romanicus"
      SPECIES["europaeus"]   = "Erinaceus_europaeus"
   }(/^Er/){
      print $1 "_0\t" SPECIES[$2] "\n" $1 "_1\t" SPECIES[$2]
   }' $POPMAP > popmap.txt
fi

if [ ! -d nexus ]; then mkdir nexus; fi

# The loci2nex.awk script selects loci with more than 2 parsimony-informative
# sites. And it excludes the hybrid individual from the alignments, because
# the multispecies model assumes reproductive isolation among species.
if ! ls nexus | grep -q ".nex"; then
   gawk -f loci2nex.awk $LOCI
   mv *.nex nexus/
fi

# To remove loci with signs of recombination, I use an existent perl script.
# by Jeffrey Ross-Ibarra. It only accepts fasta input.

touch checkpoints
if ! grep -q recombination checkpoints; then
   echo "Checking for recombination" >> checkpoints
   if [ ! -d recombinant ]; then mkdir recombinant; fi
   for i in `ls -1 nexus`; do
      gawk '(/^    Er/){print ">" $1 "\n" $2}' nexus/$i > z1.fa
      RMIN=`perl RminCutter.pl -i z1.fa -v | grep RMIN: | cut -f 2 -d " "`
      if [ $RMIN -gt 0 ]; then
         mv nexus/$i recombinant/
      fi
      rm z1*
   done
fi

# I want to calculate the GC content in each alignment.
if [ ! -e GC.txt ]; then
   gawk '(/^    Er/){
      FILE = substr(FILENAME, 7, length(FILENAME) - 10)
      GC[FILE] += gsub(/[GC]/, "x", $2)
      AT[FILE] += gsub(/[AT]/, "x", $2)
   }END{
      print "#Locus\tGC%"
      for (file in GC) {
         printf "%s\t%.4f\n", file, GC[file] / (AT[file] + GC[file])
      }
   }' nexus/*.nex > GC.txt
fi

# An now the base substitution model selection:
if [ ! -d jmodeltest ]; then mkdir jmodeltest; fi
for i in `ls -1 nexus`; do
   if [ ! -e jmodeltest/$(basename $i .nex).out ]; then
      gawk '(/dimensions/){
         split($2, NTAX,  /=/)
         split($3, NCHAR, /[=;]/)
         print NTAX[2] " " NCHAR[2]
      }(/^    Er/){
         split($1, NAME, /_/)
         printf "% -10s  %s\n", NAME[1] "_" NAME[3], $2
      }' nexus/$i > z1.phy

      java -jar $JMODELTEST_HOME/jModelTest.jar -d z1.phy -AIC -BIC -DT -dLRT -o jmodeltest/$(basename $i .nex).out -s 11 -tr 12
      rm z1.phy
   fi
done

# RevBayes does need a species map for each locus.
if [ ! -d species_maps ]; then mkdir species_maps; fi
for i in `ls -1 nexus/`; do
   if [ ! -e species_maps/$i.map ]; then
      grep $( gawk '(/^    Er/){PATTERN = PATTERN " -e " $1}END{print PATTERN}' nexus/$i ) popmap.txt > species_maps/$(basename $i .nex).map
   fi
done

# I find a limit of 8188 characters in the possible length of a command to define
# a vector of loci names. There are other very good reasons to reduce the number of loci
# analysed. Start thinking how to select the most informative sites.
