#!/bin/bash
#
#				2020-07-19
#				==========
#
# Kristyna made the file "data_description_dryad.txt" to document the files uploaded
# to dryad. She extensively documented several aspects of the samples, like their
# barcodes, the data sets they are included in, geographical coordinates, species
# assignation... There are a couple of minor details I would like to edit in that
# file before submission to dryad. Because of how important these metadata are, I
# do not dare editing it with bare hands, and I use this script instead, hoping my
# mistakes will get recorded and corrected later, if needed.
#
# This is a summary of the file:
#
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# Sample ID       Alternative Sample ID   Barcode Dataset 1       Dataset 2       Dataset 3       Dataset 4       Dataset 5       Dataset 6       Coordinates             species determination accor
# 1       Er26_JUG4       ACCAT   yes     yes     no      yes     yes     yes     41.11   20.82   Erinaceus roumanicus    Balkan lineage of E. roumanicus
# 2       Er27_SK32       ACGTA   yes     yes     no      yes     yes     yes     48.24   16.95   Erinaceus roumanicus    Balkan lineage of E. roumanicus
# 3       Er28_112        ACTGC   yes     yes     no      yes     yes     yes     50.12   14.42   Erinaceus roumanicus    Balkan lineage of E. roumanicus
# 4       Er29_122        AGTCA   yes     yes     no      yes     yes     yes     50.12   14.42   Erinaceus europaeus     Apennine lineage of E. europaeus
# 5       Er30_79 ATGCT   yes     yes     no      yes     yes     yes     54.77   -1.59   Erinaceus europaeus     Iberian lineage of E. europaeus
# 6       Er31_453        CCTTG   yes     yes     no      yes     yes     yes     53.55   IX.99   Erinaceus europaeus     Apennine lineage of E. europaeus
# ...
# 45      Er70_RMN42      GTCAC   yes     yes     no      yes     yes     yes     46.83   24.74   Erinaceus roumanicus    Balkan lineage of E. roumanicus
# 46      Er71_SAR2       GTGTG   yes     yes     no      yes     no      yes     40.31   IX.19   Erinaceus europaeus     Iberian/Apennine lineage of E. europaeus
# 47      Er72_SIE1B      GTTGT   yes     yes     no      yes     yes     yes     38.03   14.II   Erinaceus europaeus     Iberian lineage of E. europaeus
# 48      Er73_SNG1       TAATG   no      no      no      no      no      no      13.65   -16.45  Atelerix albiventris    no lineage
# 49      Er74_SP16       TACGT   yes     yes     no      yes     yes     yes     41.67   II.26   Erinaceus europaeus     Iberian lineage of E. europaeus
# 50      Er75_TRC2A      TAGCA   yes     yes     yes     no      no      yes     36.14   36.07   Erinaceus concolor      no lineage
#
# RAD single read
# Sample ID       Alternative Sample ID   Barcode Dataset 1       Dataset 2       Dataset 3       Dataset 4       Dataset 5       Dataset 6       Coordinates
# 51      Er76_Ch8f†      AAAAA   yes     no      no      yes     yes     yes     56.07   38.33   Erinaceus europaeus     Asian lineage of E. roumanicus
# 52      Er77_Ch18f†     AACCC   yes     no      no      yes     yes     yes     56.07   38.33   Erinaceus europaeus     F1 Asian lineage x Apennine lineage
# 53      Er78_Ch11       AAGGG   yes     no      yes     yes     yes     yes     56.07   38.33   Erinaceus roumanicus    Asian lineage of E. roumanicus
# 54      Er79_Ch19       AATTT   yes     no      yes     yes     yes     yes     55.61   38.54   Erinaceus roumanicus    Asian lineage of E. roumanicus
# 55      Er80_Ch23       ACACG   yes     no      yes     yes     yes     yes     56.07   38.33   Erinaceus europaeus     Apennine lineage of E. europaeus
# 56      Er81_Ch44†      ACCAT   yes     no      no      yes     yes     yes     56.07   38.33   Erinaceus roumanicus    Asian lineage of E. roumanicus
# ...
# 85      Er110_Krasngr   GATCG   no      no      yes     no      no      yes     55.83   37.33   Erinaceus europaeus     no lineage
# 86      Er111_Saratov111        GCATT   no      no      yes     no      no      yes     50.78   46.92   Erinaceus roumanicus    Asian lineage of E. roumanicus
# 87      Er112_Samara171 GCCGG   yes     no      yes     yes     yes     yes     53.43   49.67   Erinaceus roumanicus    Asian lineage of E. roumanicus
# 88      Er113_KaluE6    GCGCC   yes     no      yes     yes     yes     yes     53.63   35.75   Erinaceus roumanicus    Asian lineage of E. roumanicus
# 89      Er114_Penza171  GCTAA   yes     no      yes     yes     yes     yes     53.23   45.25   Erinaceus roumanicus    Asian lineage of E. roumanicus
# 90      Er115_VortsaE12†        GGAAG   yes     no      no      yes     yes     yes     58.11   52.15   Erinaceus europaeus     Apennine lineage of E. europaeus
#
# † hybrids.
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# The problems I see with this data set are related to how difficult it would be to
# parse it automatically:
#
# 1. There are two tables, with their headers, instead of only one.
# 2. The information of whether it is a hybrid or not is encoded with a strange
#    character, as a footnote.
# 3. Some column names are too long, and all they need some explanation.
# 4. Some coordinates are in roman numbers! I'm not going to change that...
# 5. Coordinates span two columns.

ORIGINAL=/data/kristyna/hedgehog/data_2020/data_description_dryad.txt

# The strange character is the dagger character, that can be matched with \xE2\x80\xA0.

if [ ! -e ../../data/dryad/metadata.txt ]; then
   gawk -v FS="\t" 'BEGIN{
      print "# Sample_ID1: sample identification number."
      print "# Sample_ID2: sample name used in most files and output."
      print "# Barcode: inline codeword that identifies the sample in multiplexed sequencing."
      print "# Dataset_X: whether the sample is included in data set X, with X between 1 and 6."
      print "# Coordinates: geographic coordinates of sampling location."
      print "# Mt-species: species assignment according to mitochondrial D-loop sequence."
      print "# Lineage: species and lineage assignment according to Admixture analysis."
      print "# Hybrid: whether hybrid ancestry (E. roumanicus x E. europaeus) is detected in Admixture analysis."
      print "# Sequencing: sequencing mode (paired-end or single-read).\n"
      print "Sample_ID1\tSample_ID2\tBarcode\tDataset_1\tDataset_2\tDataset_3\tDataset_4\tDataset_5\tDataset_6\tCoordinates\tMt-species\tLineage\tHybrid\tSequencing"
      SEQMODE="paired-end"
   }(/^RAD single read/){
      SEQMODE="single-read"
   }(/Er/){
      if (/\xE2\x80\xA0/) {
         HYBRID="yes"
         gsub(/\xE2\x80\xA0/,"")
      } else {
         HYBRID="no"
      }
      print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "," $11 "\t" $12 "\t" $13 "\t" HYBRID "\t" SEQMODE
   }' $ORIGINAL > ../../data/dryad/metadata.txt
fi
