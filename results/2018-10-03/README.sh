#!/bin/bash
#
#				2018-10-03
#				==========
#
# On 2018-06-30 I implemented my own script to infer local ancestry of admixed individuals
# from the output of Admixture. After Kristýna has repeated the Admixture analysis,
# with newly filtered vcf files (see 2018-09-26) I want to infer local ancestry blocks with
# more standard tools. The current programs to infer chromosomal blocks of ancestry
# are designed to work with dense and phased SNP data, mostly oriented to human
# populations. Older programs are more appropriate for the simple case of a few
# putatively admixed individuals, with low-density, unphased data. My first option
# was SABER (Tang et al. 2006. https://doi.org/10.1086/504302), which is implemented
# in R and accepts unphased, genotype data. However, I was not able to install it,
# and I am waiting for a response from the author about it. The second program is
# called LAMP (Sankararaman et al. 2008. https://doi.org/10.1016/j.ajhg.2007.09.022).
# This is the one I will try here.
#
# Each presumably admixed individual should be analysed separately. There is the hybrid
# individual, with E. roumanicus background and ~25% E. europaeus ancestry, another
# E. roumanicus sample with 4.8% E. europaeus ancestry (Er55_AU7), and the two
# additional E. roumanicus individuals with a small percentage of putative E. concolor
# ancestry. Unfortunately, chromosomes (i.e., contigs) must be analysed separately, in
# different folders. Oops! In each individual- and contig-specific folder, I should
# have the genotypes, frequencies and the positions files. It would be nice to output
# those files along the same reading of the vcf file. I realize that I don't actually
# need the Admixture output to run LAMP, since I know what individuals belong to each
# population, and I can estimate population-specific allele frequencies from the vcf
# file. I will take advantage of the binary presence flags to speed it up. The only
# thing I need from the Admixture output is the proportions of ancestries in each
# individual.
#
# After a few runs, I see that the results depende quite a bit on the amount of recombination
# expected. And we don't really know how long ago the admixture happened in the lineage
# of some individuals. I would like to run the whole analysis for several values of the
# number of generations since admixture. What I do know, or at least can estimate, is
# the average recombination rate between a pair of SNPs per generation. Er37_SK27 is a
# backcross: one of its parents was a hybrid, and it is the recombination in that parent
# what produced all the blocks of alternative ancestry. Note that it is not necessary
# to know the phase in order to count the number of crossovers in a single parent. From
# previous runs of LAMP, I approximately count 5 clear breakpoints within the first 50
# contigs, which mean a recombination rate per base per generation of 1.0e-8.

ADMIXTURE_DIR=/data/kristyna/hedgehog/results_2018/05-06-2018/r10e10c4
VCF_FILE=../2018-09-26/ErinMaxMiss63_r10e10c4.vcf
ADMIXED=(  Er37_SK27   Er55_AU7  Er26_JUG4  Er27_SK32 )
ANCEST1=( roumanicus roumanicus roumanicus roumanicus )
ANCEST2=(  europaeus  europaeus   concolor   concolor )

# This is the minimum number of SNPs per contig to analyse the contig:
MIN_LOCI=50

# And this, the number of largest contigs to represent graphically.
MAX_CONTIGS=100

# I want to runt LAMP for a different set of values of the number-of-generations
# parameter, for each individual. Because I iterate over individuals below, it
# is not straight forward to tell the script what values to iterate over for each
# individual. It can be done with "eval", but indirect references are supposed
# to be preferred. E.g.: "A=B; B=3; echo ${!A}" outputs 3. However indirect references
# do now work well with arrays, because "${!name[@]}" expands to the list of keys
# in the array 'name', not to the array it may be indirectly referring to. Maybe the
# simplest solution is to use an associative array, with the list of values for each
# individual introduced as a string, that can later be split.

declare -A GENERATIONS
GENERATIONS[Er37_SK27]='2;3;4;5;6;7;8;9;10;25;50'
GENERATIONS[Er55_AU7]='2;5;10;25;50;100;250;500;1000;2500;5000;10000;25000;50000;100000;250000;500000;1000000'
GENERATIONS[Er26_JUG4]='2;5;10;25;50;100;250;500;1000;2500;5000;10000;25000;50000;100000;250000;500000;1000000'
GENERATIONS[Er27_SK32]='2;5;10;25;50;100;250;500;1000;2500;5000;10000;25000;50000;100000;250000;500000;1000000'

if [ ! -e popmap.txt ]; then
   cp ../2018-09-26/popmap.txt ./
fi

if [ ! -e Q3.txt ]; then
   paste $ADMIXTURE_DIR/out.012.indv $ADMIXTURE_DIR/erinaceus_41_r10e10c4.3.Q > Q3.txt
fi

for ind in ${ADMIXED[@]}; do
   if [ ! -d $ind ]; then mkdir $ind; fi
done

if [ ! -e contig.dict ]; then
   gawk '(/^##contig=<ID=/){split($1,A,/[=,>]/); print A[3] "\t" A[5]}' $VCF_FILE > contig.dict
fi

for contig in $(cut -f 1 contig.dict); do
   for ind in ${ADMIXED[@]}; do
      if [ ! -d $ind/$contig ]; then mkdir $ind/$contig; fi
   done
done

# This generates the data files required by lamp.
if [ ! -e lamp_preparation.log ]; then
   gawk -v IND="Er37_SK27,Er55_AU7,Er26_JUG4,Er27_SK32" \
        -v POP1="roumanicus,roumanicus,roumanicus,roumanicus" \
        -v POP2="europaeus,europaeus,concolor,concolor" \
        -f prepare_lamp.awk popmap.txt $VCF_FILE > lamp_preparation.log
fi

# This is meant to remove from analysis contigs without data.
find . -type d -empty -delete

for ind in ${ADMIXED[@]}; do
   # This creates a joint site frequency spectrum for the two ancestral
   # populations of each admixed individual. Just out of curiosity.
   if [ ! -e $ind/joint_SFS.txt ]; then
      find $ind/ -mindepth 1 -name 'NW*' -type d -exec paste '{}'/Anc1.P.txt '{}'/Anc2.P.txt \; | \
      gawk '{
         F[sprintf("%.1f\t%.1f", $1,$2)]++
      }END{
         for (x=1.0;x>=0.0;x-=0.1) {
            S=""
            for (y=0.0;y<=1.0;y+=0.1) S = S "\t" F[sprintf("%.1f\t%.1f",x,y)] + 0
            print S
         }
      }' > $ind/joint_SFS.txt
   fi

   # I should also remove now contigs with fewer positions than the minimum.
   for contig in $(find $ind -type d -name 'NW*'); do
      if [ $(cat $contig/positions.txt | wc -l) -lt $MIN_LOCI ]; then
         rm -r $contig
      fi
   done

   GEN_STRING=${GENERATIONS[$ind]}
   GEN_ARRAY=( ${GEN_STRING//;/ } )   # This substitutes all ";" for " " in GEN_STRING
   for gen in ${GEN_ARRAY[@]}; do
#      if [ ! -e $(printf "%s/summary_gen%06i.txt" $ind $gen) ]; then
         ./run_lamp.sh $ind $gen 1> run_lamp_$ind.log 2> run_lamp_$ind.err &
#      fi
   done
   wait
done

LONGEST=$(sort -nrk 2,2 contig.dict | head -n 1 | cut -f 1)
for ind in ${ADMIXED[@]}; do
   GEN_STRING=${GENERATIONS[$ind]}
   GEN_ARRAY=( ${GEN_STRING//;/ } )
   for gen in ${GEN_ARRAY[@]}; do
      if [ ! -e $(printf "%s/mosaic%06i.png" $ind $gen) ]; then
         SCALE=$(echo "$(file $(printf "%s/%s/gen%06i/tmp/chr.png" $ind $LONGEST $gen) | cut -d ' ' -f 5) / $(sort -nrk 2,2 contig.dict | head -n 1 | cut -f 2)" | bc -l)
         HEIGHT=$(file $(printf "%s/%s/gen%06i/tmp/chr.png" $ind $LONGEST $gen) | cut -d ' ' -f 7)
         if [ ! -e $(printf "%s/args%06i.txt" $ind $gen) ]; then
            for contig in $(sort -nrk 2,2 contig.dict | head -n $MAX_CONTIGS | cut -f 1); do
               if [[ -e $(printf "%s/%s/gen%06i/tmp/chr.png" $ind $contig $gen) ]] && [[ ! -e $(printf "%s/%s/gen%06i/tmp/chr_scaled.png" $ind $contig $gen) ]]; then
                  WIDTH=$(echo -e "scale=0\n$(grep $contig contig.dict | cut -f 2) * $SCALE" | bc -l)
                  convert $(printf "%s/%s/gen%06i/tmp/chr.png" $ind $contig $gen) -resize ${WIDTH}x${HEIGHT%,}! -quality 100 $(printf "%s/%s/gen%06i/tmp/chr_scaled.png" $ind $contig $gen)
                  echo $(printf "%s/%s/gen%06i/tmp/chr_scaled.png" $ind $contig $gen) >> $(printf "%s/args%06i.txt" $ind $gen)
                  echo $contig >> z1
               fi
            done
            printf -- "-append %s/mosaic%06i.png" $ind $gen >> $(printf "%s/args%06i.txt" $ind $gen)
         fi
         xargs --arg-file=$(printf "%s/args%06i.txt" $ind $gen) convert -background black
      fi
#      if [ -e $(printf "%s/args%06i.txt" $ind $gen) ]; then rm $(printf "%s/args%06i.txt" $ind $gen); fi
   done
done

if [[ ! -e summary.txt ]] || [[ -n $(find . -name 'summary_gen*' -newer summary.txt) ]]; then
   echo -e "#Sample\tAncestral_1\tAncestral_2\tPercent1_(SNPs)\tPercent2_(SNPs)\tPercent1_(bases)\tPercent2_(bases)\tGenerations\tLogLikelihood" > summary.txt
   for i in 0 1 2 3; do
      GEN_STRING=${GENERATIONS[${ADMIXED[$i]}]}
      GEN_ARRAY=( ${GEN_STRING//;/ } )
      for gen in ${GEN_ARRAY[@]}; do
         LOGL=`find $(printf "%s/*/gen%06i" ${ADMIXED[$i]} $gen) -name lamp.log -exec tail -n 1 '{}' \; | gawk '{split($4,A,/=/); S += A[2]}END{printf("%.4f\n", S)}'`
         gawk -v IND=${ADMIXED[$i]} -v ANC1=${ANCEST1[$i]} -v ANC2=${ANCEST2[$i]} -v GEN=$gen -v LOGL=$LOGL '{
            WSUM_CONTIG_1 += $2 * $4
            WSUM_CONTIG_2 += $2 * $5
            TOTAL_CONTIG  += $2
            WSUM_SNP_1    += $3 * $4
            WSUM_SNP_2    += $3 * $5
            TOTAL_SNP     += $3
         }END{
            printf("%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t% 8i\t% 13.4f\n", IND, ANC1, ANC2, WSUM_SNP_1 / TOTAL_SNP, WSUM_SNP_2 / TOTAL_SNP, \
              WSUM_CONTIG_1 / TOTAL_CONTIG, WSUM_CONTIG_2 / TOTAL_CONTIG, GEN, LOGL)
         }' $(printf "%s/summary_gen%06i.txt" ${ADMIXED[$i]} $gen) >> summary.txt
      done
   done
fi
