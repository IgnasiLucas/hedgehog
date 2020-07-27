#!/bin/bash
#
#				2020-07-22
#				==========
#
# In 2020-07-14 I prepared the tar files to submit to Dryad. I took care
# to keep the files in an organized way, so that each of the three data sets
# can be downloaded separately, and create separate folders. The three
# tar files ended up being 27, 44 and 46 GB. Because I'm working remotely,
# and I don't have illimited internet connexion, I requested from Dryad
# an ftp upload option. Even though I had told them how big the files are,
# they told me I can use ftp, but only to upload files up to 10GB in size.
#

DATADIR1='/data/kristyna/hedgehog/data_2016/fastq'
DATADIR2='/data/kristyna/hedgehog/data_2018/fastq'
DATADIR3='/data/kristyna/hedgehog/data_2020'

if [ ! -d ../../data/dryad ]; then mkdir ../../data/dryad; fi

EUR=(Er26_JUG4 Er27_SK32 Er28_112   Er29_122  Er30_79    Er31_453  Er32_183   Er33_211
     Er34_197  Er35_209  Er36_SK24  Er37_SK27 Er38_LB1   Er39_PL1  Er40_M2    Er41_GR36
     Er42_GR35 Er43_SL7  Er44_VOJ1  Er45_BLG3 Er46_RMN7  Er47_CR4  Er48_BH16  Er49_GR38
     Er50_R3   Er51_436  Er52_451   Er53_ASR7 Er54_AU1   Er55_AU7  Er56_AZ5   Er57_COR4
     Er58_FI7  Er59_FR1  Er60_GR5   Er61_GR87 Er62_GR95  Er63_IR6  Er64_IS1   Er65_IS25
     Er66_IT3  Er67_IT5  Er68_PRT1B Er69_R2   Er70_RMN42 Er71_SAR2 Er72_SIE1B Er73_SNG1
     Er74_SP16 Er75_TRC2A)

RUS=(Er100_ZBS133    Er101_ZBS135 Er102_ZBS1675   Er103_ZBSTvE8   Er104_Kal121  Er105_Vor151
     Er106_NGE1      Er107_Kos27  Er108_VLD2013   Er109_NVB       Er110_Krasngr Er111_Saratov111
     Er112_Samara171 Er113_KaluE6 Er114_Penza171  Er115_VortsaE12 Er76_Ch8f     Er77_Ch18f
     Er78_Ch11       Er79_Ch19    Er80_Ch23       Er81_Ch44       Er82_Ch53     Er83_Ea108
     Er84_Ea091      Er85_TuE2    Er86_Tv1302     Er87_Tv161      Er88_PTZ092   Er89_PTZ1301
     Er90_PTZ1303    Er91_GR24    Er92_GR28       Er93_ZBS073     Er94_ZBS093   Er95_ZBS104
     Er96_ZBS1011    Er97_ZBS1015 Er98_ZBS103.111 Er99_ZBS122)

# ========================================================================================
#                                     2016 data set
# ========================================================================================

if ! $(ls -l ../../data/dryad/ | grep -q fastq2016); then
   TAR=1
   touch ../../data/dryad/fastq2016_$TAR.tar
   for i in ${EUR[@]}; do
      TAR_SIZE=$(ls -l ../../data/dryad/fastq2016_$TAR.tar | cut -d " " -f 5)
      FQ1_SIZE=$(ls -l $DATADIR1/${i}_1.fastq.gz | cut -d " " -f 5)
      FQ2_SIZE=$(ls -l $DATADIR1/${i}_2.fastq.gz | cut -d " " -f 5)
      if [ $(( $TAR_SIZE + $FQ1_SIZE + $FQ2_SIZE )) -lt 10000000000 ]; then
         tar -rf ../../data/dryad/fastq2016_$TAR.tar -P --transform "s,$DATADIR1,fastq2016," $DATADIR1/${i}_1.fastq.gz $DATADIR1/${i}_2.fastq.gz
      else
         TAR=$(( TAR + 1 ))
         tar -cf ../../data/dryad/fastq2016_$TAR.tar -P --transform "s,$DATADIR1,fastq2016," $DATADIR1/${i}_1.fastq.gz $DATADIR1/${i}_2.fastq.gz
      fi
   done
   if ! $(tar -tf ../../data/dryad/fastq2016_$TAR.tar | grep -q md5sum); then
      if [ ! -e fastq2016_md5sum.txt ]; then
         md5sum $DATADIR1/Er*.fastq.gz | sed "s,$DATADIR1,fastq2016," > fastq2016_md5sum.txt
      fi
      tar -rf ../../data/dryad/fastq2016_$TAR.tar fastq2016_md5sum.txt
   fi
fi

# ========================================================================================
#                                     2018 data set
# ========================================================================================

if ! $(ls -l ../../data/dryad/ | grep -q fastq2018); then
   TAR=1
   touch ../../data/dryad/fastq2018_$TAR.tar
   for i in ${EUR[@]}; do
      TAR_SIZE=$(ls -l ../../data/dryad/fastq2018_$TAR.tar | cut -d " " -f 5)
      FQ1_SIZE=$(ls -l $DATADIR2/${i}_1.fastq.gz | cut -d " " -f 5)
      FQ2_SIZE=$(ls -l $DATADIR2/${i}_2.fastq.gz | cut -d " " -f 5)
      if [ $(( $TAR_SIZE + $FQ1_SIZE + $FQ2_SIZE )) -lt 10000000000 ]; then
         tar -rf ../../data/dryad/fastq2018_$TAR.tar -P --transform "s,$DATADIR2,fastq2018," $DATADIR2/${i}_1.fastq.gz $DATADIR2/${i}_2.fastq.gz
      else
         TAR=$(( TAR + 1 ))
         tar -cf ../../data/dryad/fastq2018_$TAR.tar -P --transform "s,$DATADIR2,fastq2018," $DATADIR2/${i}_1.fastq.gz $DATADIR2/${i}_2.fastq.gz
      fi
   done
   if ! $(tar -tf ../../data/dryad/fastq2018_$TAR.tar | grep -q md5sum); then
      if [ ! -e fastq2018_md5sum.txt ]; then
         md5sum $DATADIR2/Er*.fastq.gz | sed "s,$DATADIR2,fastq2018," > fastq2018_md5sum.txt
      fi
      tar -rf ../../data/dryad/fastq2018_$TAR.tar fastq2018_md5sum.txt
   fi
fi

# ========================================================================================
#                                     2020 data set
# ========================================================================================

if ! $(ls -l ../../data/dryad/ | grep -q fastq2020); then
   TAR=1
   touch ../../data/dryad/fastq2020_$TAR.tar
   for i in ${RUS[@]}; do
      TAR_SIZE=$(ls -l ../../data/dryad/fastq2020_$TAR.tar | cut -d " " -f 5)
      FQ1_SIZE=$(ls -l $DATADIR3/${i}.fastq.gz | cut -d " " -f 5)
      if [ $(( $TAR_SIZE + $FQ1_SIZE )) -lt 10000000000 ]; then
         tar -rf ../../data/dryad/fastq2020_$TAR.tar -P --transform "s,$DATADIR3,fastq2020," $DATADIR3/${i}.fastq.gz
      else
         TAR=$(( TAR + 1 ))
         tar -cf ../../data/dryad/fastq2020_$TAR.tar -P --transform "s,$DATADIR3,fastq2020," $DATADIR3/${i}.fastq.gz
      fi
   done
   if ! $(tar -tf ../../data/dryad/fastq2020_$TAR.tar | grep -q md5sum); then
      if [ ! -e fastq2020_md5sum.txt ]; then
         md5sum $DATADIR3/Er*.fastq.gz | sed "s,$DATADIR3,fastq2020," > fastq2020_md5sum.txt
      fi
      tar -rf ../../data/dryad/fastq2020_$TAR.tar fastq2020_md5sum.txt
   fi
fi
