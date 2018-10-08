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

ADMIXTURE_DIR=/data/kristyna/hedgehog/results_2018/05-06-2018/r10e10c4
VCF_FILE=../2018-09-26/ErinMaxMiss63_r10e10c4.vcf
ADMIXED=(  Er37_SK27   Er55_AU7  Er26_JUG4  Er27_SK32 )
ANCEST1=( roumanicus roumanicus roumanicus roumanicus )
ANCEST2=(  europaeus  europaeus   concolor   concolor )
# This is the minimum number of SNPs per contig to analyse the contig.
MIN_LOCI=50
#GENERATIONS=( 2 10 20000 20000 )
#GENERATIONS=( 3 100 10000 10000 )
GENERATIONS=( 4 200 5000 5000 )

if [ ! -e popmap.txt ]; then
   cp ../2018-09-26/popmap.txt ./
fi

if [ ! -e Q3.txt ]; then
   paste $ADMIXTURE_DIR/out.012.indv $ADMIXTURE_DIR/erinaceus_41_r10e10c4.3.Q > Q3.txt
fi

for ind in Er37_SK27 Er55_AU7 Er26_JUG4 Er27_SK32; do
   if [ ! -d $ind ]; then mkdir $ind; fi
done

for contig in `grep '^##contig=<ID=' $VCF_FILE | sed -r 's/(##contig=<ID=|,length=[0-9]+>)//g'`; do
   for ind in Er37_SK27 Er55_AU7 Er26_JUG4 Er27_SK32; do
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

find . -type d -empty -delete

if [ ! -e contig.dict ]; then
   gawk '(/^##contig=<ID=/){split($1,A,/[=,>]/); print A[3] "\t" A[5]}' $VCF_FILE > contig.dict
fi

for ind in Er37_SK27 Er55_AU7 Er26_JUG4 Er27_SK32; do
   # This creates a joint site frequency spectrum for the two ancestral
   # populations of each admixed individual. Just out of curiosity.
   if [ ! -e $ind/joint_SFS.txt ]; then
      find $ind/ -mindepth 1 -type d -exec paste '{}'/Anc1.P.txt '{}'/Anc2.P.txt \; | \
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

   if [ ! -e $ind/configuration.txt ]; then
      echo "populations=2"                    > $ind/configuration.txt
      echo "genofile=geno.txt"               >> $ind/configuration.txt
      echo "posfile=positions.txt"           >> $ind/configuration.txt
      echo "pfile=Anc1.P.txt,Anc2.P.txt"     >> $ind/configuration.txt
      ASZ1=`grep -E "^ ?$ind" lamp_preparation.log | gawk '{printf("%.0f",$5)}'`
      ASZ2=`grep -E "^ ?$ind" lamp_preparation.log | gawk '{printf("%.0f",$6)}'`
      echo "ancestralsamplesize=$ASZ1,$ASZ2" >> $ind/configuration.txt
      ALPH=`grep "^$ind" Q3.txt | gawk '{print $3}'`
      BETA=0`echo "1.0 - $ALPH" | bc -l`
      echo "alpha=$ALPH,$BETA"               >> $ind/configuration.txt
      GEN=`for i in 0 1 2 3; do echo -e "${ADMIXED[$i]}\t${GENERATIONS[$i]}"; done | grep $ind | cut -f 2`
      echo "generations=$GEN"                >> $ind/configuration.txt
      echo "recombrate=1.0e-9"               >> $ind/configuration.txt
   fi
   for contig in `find $ind/ -mindepth 1 -type d`; do
      cd $contig
         cp ../configuration.txt ./configuration.txt
         if [ `cat positions.txt | wc -l` -ge $MIN_LOCI ] && [ ! -e ancestry_lamp4.out ]; then
            lamp configuration.txt 1> lamp.log 2> lamp.err
         fi
      cd ../..
   done
   if [ ! -e $ind/summary.txt ]; then
      ./summarize_lamp.sh $ind contig.dict | sort -nrk 2,2 > $ind/summary.txt
   fi
done

if [ ! -e summary.txt ]; then
   echo -e "#Sample\tAncestral_1\tAncestral_2\tPercent1_(SNPs)\tPercent2_(SNPs)\tPercent1_(bases)\tPercent2_(bases)" > summary.txt
fi
echo "" >> summary.txt
for i in 0 1 2 3; do
   LOGL=`find ${ADMIXED[$i]} -name lamp.log -exec tail -n 1 '{}' \; | gawk '{split($4,A,/=/); S += A[2]}END{printf("%.4f\n", S)}'`
   gawk -v IND=${ADMIXED[$i]} -v ANC1=${ANCEST1[$i]} -v ANC2=${ANCEST2[$i]} -v GEN=${GENERATIONS[$i]} -v LOGL=$LOGL '{
      WSUM_CONTIG_1 += $2 * $4
      WSUM_CONTIG_2 += $2 * $5
      TOTAL_CONTIG  += $2
      WSUM_SNP_1    += $3 * $4
      WSUM_SNP_2    += $3 * $5
      TOTAL_SNP     += $3
   }END{
      printf("%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t% 6i\t% 8.4f\n", IND, ANC1, ANC2, WSUM_SNP_1 / TOTAL_SNP, WSUM_SNP_2 / TOTAL_SNP, \
        WSUM_CONTIG_1 / TOTAL_CONTIG, WSUM_CONTIG_2 / TOTAL_CONTIG, GEN, LOGL)
   }' ${ADMIXED[$i]}/summary.txt >> summary.txt
done
