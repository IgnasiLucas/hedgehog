#!/bin/bash
#
#  Usage: run_lamp.sh <ind> <gen> <min_loci>
#

WORKDIR=/data/joiglu/hedgehog/results/2018-10-03

if [ ! -e $(printf "%s/configuration_%06i.txt" $1 $2) ]; then
   echo "populations=2"                          > $(printf "%s/configuration_%06i.txt" $1 $2)
   echo "genofile=../geno.txt"                  >> $(printf "%s/configuration_%06i.txt" $1 $2)
   echo "posfile=../positions.txt"              >> $(printf "%s/configuration_%06i.txt" $1 $2)
   echo "pfile=../Anc1.P.txt,../Anc2.P.txt"     >> $(printf "%s/configuration_%06i.txt" $1 $2)
   ASZ1=`grep -E "^ ?$1" lamp_preparation.log | gawk '{printf("%.0f",$5)}'`
   ASZ2=`grep -E "^ ?$1" lamp_preparation.log | gawk '{printf("%.0f",$6)}'`
   echo "ancestralsamplesize=$ASZ1,$ASZ2"       >> $(printf "%s/configuration_%06i.txt" $1 $2)
   ALPH=`grep "^$1" Q3.txt | gawk '{print $3}'`
   BETA=0`echo "1.0 - $ALPH" | bc -l`
   echo "alpha=$ALPH,$BETA"                     >> $(printf "%s/configuration_%06i.txt" $1 $2)
   echo "generations=$2"                        >> $(printf "%s/configuration_%06i.txt" $1 $2)
   echo "recombrate=1.0e-8"                     >> $(printf "%s/configuration_%06i.txt" $1 $2)
fi

for contig in $(find $1/ -mindepth 1 -name 'NW*' -type d); do
#   if [ -d $contig ]; then
      cd $contig
         if [ ! -d $(printf "gen%06i" $2) ]; then mkdir $(printf "gen%06i" $2); fi
         cd $(printf "gen%06i" $2)
            if [ ! -e $(printf "configuration_%06i.txt" $2) ]; then
               cp $(printf "../../configuration_%06i.txt" $2) $(printf "configuration_%06i.txt" $2)
            fi
            if [ ! -e ancestry_lamp4.out ]; then
               lamp $(printf "configuration_%06i.txt" $2) 1> lamp.log 2> lamp.err
               generategraph.sh ancestry_lamp4.out
            fi
      cd $WORKDIR
#   else
#      echo "$contig does not seem to exist in "$(pwd)
#   fi
done

if [[ ! -e $(printf "%s/summary_gen%06i.txt" $1 $2) ]] || \
   [[ -n $(find $(printf "%s/gen%06i" $1 $2) -name ancestry_lamp4.out -newer $(printf "%s/summary_gen%06i.txt" $1 $2)) ]]; then
   ./summarize_lamp.sh $1 $2 contig.dict | sort -nrk 2,2 > $(printf "%s/summary_gen%06i.txt" $1 $2)
fi

