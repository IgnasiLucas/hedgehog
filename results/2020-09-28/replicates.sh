#!/bin/bash

P2=$1
P3=$2
FIRST=$3
LAST=$4
OUTDIR=$(printf "%s_%s" $P2 $P3)
if [ $P2 == "balkan" ];    then POP2="E_roumanicus_Balkan"; fi
if [ $P2 == "asian"  ];    then POP2="E_roumanicus_Asian"; fi
if [ $P3 == "iberian" ];   then POP3="E_europaeus_Iberian"; fi
if [ $P3 == "apennine" ];  then POP3="E_europaeus_Apennine"; fi
if [ $P3 == "amurensis" ]; then POP3="E_amurensis"; fi
for i in $(seq $FIRST $LAST); do
   POPMAP=$(printf "%s/popmap%05u.txt" $OUTDIR $i)
   if [ ! -e $(printf "%s/rep%05u.freq" $OUTDIR $i) ]; then
      gawk --bignum \
        -v P1="E_concolor" \
        -v P2=$POP2 \
        -v P3=$POP3 \
        -v OUTGROUP="Hemiechinus" \
        -v MIN1=1 -v MIN2=1 -v MIN3=1 -v MINOUT=1 \
        -f ../../bin/freq2.awk $POPMAP EH.vcf | \
        grep -v -F 999.9999 > $(printf "%s/rep%05u.freq" $OUTDIR $i)
   fi
done
