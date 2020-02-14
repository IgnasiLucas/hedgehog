#!/bin/bash
#
#				2020-02-12
#				----------
#
# A new dataset is available. Kristýna has already mapped the reads and assembled the
# vcf file. I want to check a few things. First, if the mapping would have been better
# with different settings. And second, if the filtering of variable sites is correct.

FASTQDIR='/data/kristyna/hedgehog/data_2020'
REFERENCE='/data/kristyna/hedgehog/data_2016/reference/reference_erinaceus'
SAMPLE=(Er76_Ch8f     Er77_Ch18f   Er78_Ch11     Er79_Ch19        Er80_Ch23       Er81_Ch44    Er82_Ch53       Er83_Ea108
        Er84_Ea091    Er85_TuE2    Er86_Tv1302   Er87_Tv161       Er88_PTZ092     Er89_PTZ1301 Er90_PTZ1303    Er91_GR24
        Er92_GR28     Er93_ZBS073  Er94_ZBS093   Er95_ZBS104      Er96_ZBS1011    Er97_ZBS1015 Er98_ZBS103.111 Er99_ZBS122
        Er100_ZBS133  Er101_ZBS135 Er102_ZBS1675 Er103_ZBSTvE8    Er104_Kal121    Er105_Vor151 Er106_NGE1      Er107_Kos27
        Er108_VLD2013 Er109_NVB    Er110_Krasngr Er111_Saratov111 Er112_Samara171 Er113_KaluE6 Er114_Penza171  Er115_VortsaE12)
BOWTIE2STATS='/data/kristyna/hedgehog/results_2020/2020-01-30/bowtie2_sumstat'
POPMAP='../../data/populations_2020.txt'
VCF1='/data/kristyna/hedgehog/results_2020/2020-02-02/russia.vcf'
VCF2='/data/kristyna/hedgehog/results_2020/2020-02-02/2020-02-02.recode.vcf'


# Mapping
# -------
#
# This is just formatting the kristyna's bowtie2_sumstat file.
if [ ! -e original_mapping_statistics.txt ]; then
   gawk '((FILENAME ~ /pop/) && (/^Er/)){
      POP[$1] = $2
   }((FILENAME ~ /pop/) && (/^#/)){
      SPECIES[$2] = $5
   }((FILENAME ~ /bowtie2_sumstat/) && (/^==>/)){
      sub(/_bowtie2.log/,"", $2)
      SAMPLE = $2
   }((FILENAME ~ /bowtie2_sumstat/) && (/reads; of these:/)){
      TOTAL[SAMPLE] = $1
   }((FILENAME ~ /bowtie2_sumstat/) && (/aligned 0 times/)){
      ZERO[SAMPLE] = $1
   }((FILENAME ~ /bowtie2_sumstat/) && (/aligned exactly 1 time/)){
      ONCE[SAMPLE] = $1
   }((FILENAME ~ /bowtie2_sumstat/) && (/aligned >1 times/)){
      MORE[SAMPLE] = $1
   }END{
      print "Sample\tUnmapped\tMapped_once\tAmbiguous_mapping\tTotal\tPopulation"
      for (sample in TOTAL) {
         printf("%s\t%u\t%u\t%u\t%u\tE. %s\n", sample, ZERO[sample], ONCE[sample], MORE[sample], TOTAL[sample], SPECIES[POP[sample]])
      }
   }' $POPMAP $BOWTIE2STATS > original_mapping_statistics.txt
fi

# There are four different species. But all individuals have a similar proportion
# of unmapped or multiply mapped reads. I would have expected individuals of the
# same species as the reference genome to have a higher mapping success. I wonder
# if the '--very-sensitive' setting in bowtie2 is over-mapping too many reads.
# I will repeat the mapping with a just '--sensitive' setting, to see what happens.

if [ ! -e new_mapping_statistics.txt ]; then
   for i in `seq 0 39`; do
      if [ ! -d unmapped ]; then mkdir unmapped; fi
      if [ ! -d sam ]; then mkdir sam; fi
      if [ ! -d logs ]; then mkdir logs; fi
      if [ ! -e sam/${SAMPLE[$i]}.sam ]; then
         bowtie2 --sensitive \
                 -x $REFERENCE \
                 -U $FASTQDIR/${SAMPLE[$i]}.fastq.gz \
                 --un-gz unmapped/${SAMPLE[$i]}.unmapped.fastq \
                 --rg-id ${SAMPLE[$i]} \
                 --rg "SM:"${SAMPLE[$i]} \
                 --rg "PL:ILLUMINA" \
                 --rg "DT:2020" \
                 -S sam/${SAMPLE[$i]}.sam 2> logs/${SAMPLE[$i]}.log &
      fi
   done
   wait

   gawk '((FILENAME ~ /pop/) && (/^Er/)){
      POP[$1] = $2
   }((FILENAME ~ /pop/) && (/^#/)){
      SPECIES[$2] = $5
   }((FILENAME ~ /log$/) && (/reads; of these:/)){
      SAMPLE = FILENAME
      gsub(/\.?logs?\/?/,"",SAMPLE)
      TOTAL[SAMPLE] = $1
   }((FILENAME ~ /log$/) && (/aligned 0 times/)){
      ZERO[SAMPLE] = $1
   }((FILENAME ~ /log$/) && (/aligned exactly 1 time/)){
      ONCE[SAMPLE] = $1
   }((FILENAME ~ /log$/) && (/aligned >1 times/)){
      MORE[SAMPLE] = $1
   }END{
      print "Sample\tUnmapped\tMapped_once\tAmbiguous_mapping\tTotal\tPopulation"
      for (sample in POP) {
         printf("%s\t%u\t%u\t%u\t%u\tE. %s\n", sample, ZERO[sample], ONCE[sample], MORE[sample], TOTAL[sample], SPECIES[POP[sample]])
      }
   }' $POPMAP logs/Er*.log > new_mapping_statistics.txt
fi

if [ ! -e assessment.html ]; then
   R --no-save -q -e "render_report <- function(mapping1, mapping2, vcf1, vcf2){ \
                        rmarkdown::render('assessment.Rmd', \
                           params = list(MAPPING1 = mapping1, MAPPING2 = mapping2, VCF1 = vcf1, VCF2 = vcf2), \
                           output_file = 'assessment.html')}" \
                  -e "render_report('original_mapping_statistics.txt', 'new_mapping_statistics.txt', '$VCF1', '$VCF2')"
fi
