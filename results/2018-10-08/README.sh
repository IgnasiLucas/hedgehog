#!/bin/bash
#
#				2018-10-08
#				==========
#
# This is just to update the R script that represents the Admixture's output,
# because the ordering of the individuals was not optimum. The original script
# was using an order based on the poportions of ancestry groups inferred with
# the largest k, namely k=9. However, I doubt any k above 4 makes much biological
# sense. Using ancestry proportions from only k=4 does not work either, because
# of one admixed roumanicus which gets ordered together with the western europaeus.
# The solution I came up with is to use first k=3 and then k=4 ancestry proportions
# to get them in this order: western europaeus, eastern europaeus, roumanicus and
# concolor. An alternative, suggested earlier by Barbora, is to use geographic
# coordinates. I may implement that later, when I find the coordinates.

ADMIXTURE_DIR=/data/kristyna/hedgehog/results_2018/05-06-2018/r10e10c4

for k in 2 3 4 5 6 7 8 9; do
   if [ ! -e erinaceus.$k.Q ]; then
      ln -s $ADMIXTURE_DIR/erinaceus_41_r10e10c4.$k.Q ./erinaceus.$k.Q
   fi
done

if [ ! -e out.012.indv ]; then
   ln -s $ADMIXTURE_DIR/out.012.indv ./out.012.indv
fi

if [ ! -e Rplots.pdf ]; then
   R --save < admixture_erinaceus.R
fi
