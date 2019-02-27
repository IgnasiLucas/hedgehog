#!/bin/bash
#
#				2019-01-29
#				==========
#
# In the design of the next sequencing experiment, it is important to know whether
# mapping 150 bp single reads can produce results similar to those from mapping paired
# 100 bp-long reads. The question is how much more unique are two 100-bp reads than one
# 150 bp read. If k-mer analsysis programs were ready to use 200-mers, it could be
# answered without having to map reads. Actually, I see online that people use KAT
# (Mapleson et al. 2017, Bioinformatics 33(4):574-576) with very large k-mers.
#
# Paired 100-bp reads should be just as unique as single 200-bp reads, more or less.
# I just need to count 150 or 200-bp words that are either unique or not.

REFERENCE='../../data/reference.fa'

for i in 100 150 200; do
   if [ ! -e ${i}mers ]; then
      kat hist --output_prefix ${i}mers \
               --threads 6 \
               --low 1 \
               --high 10 \
               --mer_len $i \
               --hash_size 4000000000 \
               --dump_hash \
               --verbose \
               $REFERENCE
   fi
done

if [ ! -e summary.txt ]; then
   gawk '(/^[^#]/){
      if ($1 == 1) UNIQUE[FILENAME] = $2
      if ($1 >  1) MULTIPLE[FILENAME] += $2
   }END{
      print "Length\t        Unique\t               Multiple\t            Total"
      for (k=100; k<=200; k+=50) {
         TOTAL = UNIQUE[k "mers"] + MULTIPLE[k "mers"]
         printf("%u\t% 12u (%3.2f%%)\t% 12u (%3.2f%%)\t% 12u\n", k, UNIQUE[k "mers"], 100 * UNIQUE[k "mers"] / TOTAL, MULTIPLE[k "mers"], 100 * MULTIPLE[k "mers"] / TOTAL, TOTAL)
      }
   }' *mers > summary.txt
fi

# ------------------------------------------------------------------
# Length        Unique	               Multiple	            Total
# ------------------------------------------------------------------
# 100	  2267026307 (99.29%)	    16141344 (0.71%)	  2283167651
# 150	  2271921385 (99.46%)	    12414941 (0.54%)	  2284336326
# 200	  2264956819 (99.51%)	    11182594 (0.49%)	  2276139413
# ------------------------------------------------------------------
#
# The results show that increasing the length of the reads from 100 to 200 would
# only improve the uniqueness (or mappability) of those reads from 99.29 to 99.51%.
# Not worth from my point of view.
