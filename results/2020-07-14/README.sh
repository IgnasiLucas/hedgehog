#!/bin/bash
#
#				2020-07-14
#				==========
#
# We need to prepare the fastq files and the definitive popmap files
# (anything else?) to submit to Dryad. There are three sets of fastq
# files, corresponding to the three sequencing runs. The two first sets
# are paired-end files from the same set of samples. Thus, therer are
# two files per sample in each set, and the names coincide between sets.
# For example, sample Er26_JUG4 has files Er26_JUG4_1.fastq.gz and
# Er26_JUG4_2.fastq.gz in the 2016 data set, but also in the 2018 dataset.
#
# Thus, I need to archive the three data sets separately. These are the
# paths of the fastq files:

DATADIR1='/data/kristyna/hedgehog/data_2016/fastq'
DATADIR2='/data/kristyna/hedgehog/data_2018/fastq'
DATADIR3='/data/kristyna/hedgehog/data_2020'

# They may include additional files, that I do not want to archive.
#
# I can create a different archive for each data set easily. However,
# anybody extracting the files from the second data set's archive
# could accidentaly overwrite files from the first data set, because
# of their common names. Thus, I should include a parent directory
# in the archives that allows to extract the two data sets separately.
#
# However, I do not want to include the common 'fastq' directory in
# the directory tree. In addition, to archive only the fastq.gz files
# when specifying a whole directory I would have to use the --exclude
# option for those files I don't want to include in the archive.
#
# So, what I want is for anybody extracting the content of 'fastq2016.tar.gz'
# to get a 'fastq2016' directory with the first data set, and so on.
# It would be a waste of time and space to re-create the original directory
# structure in order to create those archives. After searching for a
# solution, I found the --transform option of tar, which allows me to
# change the paths of the files that I archive using 'sed' syntax. Note
# the use of the -P option to preserve the root slash in the paths.

if [ ! -d ../../data/dryad ]; then mkdir ../../data/dryad; fi

if [ ! -e ../../data/dryad/fastq2016.tar ]; then
   tar -cf ../../data/dryad/fastq2016.tar -P --transform "s,$DATADIR1,fastq2016," $DATADIR1/Er*.fastq.gz
   md5sum $DATADIR1/Er*.fastq.gz | sed "s,$DATADIR1,fastq2016," > fastq2016_md5sum.txt
   tar -rf ../../data/dryad/fastq2016.tar fastq2016_md5sum.txt
fi

if [ ! -e ../../data/dryad/fastq2018.tar ]; then
   tar -cf ../../data/dryad/fastq2018.tar -P --transform "s,$DATADIR2,fastq2018," $DATADIR2/Er*.fastq.gz
   md5sum $DATADIR2/Er*.fastq.gz | sed "s,$DATADIR2,fastq2018," > fastq2018_md5sum.txt
   tar -rf ../../data/dryad/fastq2018.tar fastq2018_md5sum.txt
fi

if [ ! -e ../../data/dryad/fastq2020.tar ]; then
   tar -cf ../../data/dryad/fastq2020.tar -P --transform "s,$DATADIR3,fastq2020," $DATADIR3/Er*.fastq.gz
   md5sum $DATADIR3/Er*.fastq.gz | sed "s,$DATADIR3,fastq2020," > fastq2020_md5sum.txt
   tar -rf ../../data/dryad/fastq2020.tar fastq2020_md5sum.txt
fi

# I was planning to use the -C option together with wildcard '*' to select
# the files to include in the archive. I learned that tar does not use
# wildcards for creating an archive (only to list and extract) and it is
# the shell who expands the wildcards instead. But, if I did anything like
#
#   tar -czf z1.tar.gz -C $DATADIR1 *.fastq.gz
#
# By the time tar gets to the desired directory (option -C), the shell had
# already expanded *.fastq.gz in the original directory, probably resulting
# in nothing that can be found where tar is looking.
