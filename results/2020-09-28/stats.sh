#!/bin/bash

for DIR in asian_amurensis balkan_amurensis asian_iberian balkan_iberian asian_apennine balkan_apennine; do
   if [ ! -e $DIR/stats.tsv ]; then
      echo -e "#D\tf_hom\tf_d\tf_dM" > $DIR/stats.tsv
      R --no-save --slave < fabbababa.R --args $DIR/complete.freq >> $DIR/stats.tsv
   fi
   if [ ! -e $DIR/jackknife.tsv ]; then
      echo -e "#D\tf_hom\tf_d\tf_dM" > $DIR/jackknife.tsv
      R --no-save --slave < fabbababa.R --args $DIR >> $DIR/jackknife.tsv
   fi
done
