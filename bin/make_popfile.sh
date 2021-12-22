#!/bin/bash

POPFILE1='/gata/joiglu/hedgehog/data/populations.txt'
POPFILE2='/gata/joiglu/hedgehog/data/populations_2020.txt'

if [ ! -e popfile.txt ]; then
   gawk 'BEGIN{
      POP["1"] = "roumanicus"
      POP["2"] = "europaeus"
      POP["3"] = "concolor"
      POP["4"] = "roumanicus"
      POP["5"] = "hemiechinus"
      POP["6"] = "atelerix"
   }{print $1 "\t" POP[$2]}' $POPFILE1 > popfile.txt

   gawk 'BEGIN{
      POP["1"] = "europaeus"
      POP["2"] = "roumanicus"
      POP["3"] = "concolor"
      POP["4"] = "amurensis"
   }(/^Er/){
      print $1 "\t" POP[$2]
   }' $POPFILE2 >> popfile.txt
fi
