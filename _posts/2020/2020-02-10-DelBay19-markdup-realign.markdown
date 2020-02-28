---
layout: post
title: DelBay19 Bamfile Markdup, and Realign
date: '2020-02-10 12:04'
tags:
  - mapping
  - MAPQ
  - bowtie2
  - trimming
  - markdup
  - realign
categories:
  - WGS data analysis
---

This step will mark duplicate and realgin reads in the bam file. 

```shell

#!/bin/bash
#ACC is xxxxx.list base filename of bam filenames

REFFASTA=/home/hz269/workdir/hz269/genome/cv30.fa
PICARDDIR=/programs/bin/picard-tools
GATKDIR=/programs/GenomeAnalysisTK-3.7
TMP=/home/hz269/workdir/hz269/tmp

CPU=20

ACC=$1

# If we only want to run realignment at knownn set of sites,
# run the following only once instead of generating intervals for each BAM file

#java -Xmx1g -jar /path/to/GenomeAnalysisTK.jar \
#  -T RealignerTargetCreator \
#  -R /path/to/reference.fasta \
#  -o /path/to/output.intervals \
#  -known /path/to/indel_calls.vcf

echo Creating target intervals started
date

# Note the way Java temp files are redirected (from default /tmp directory) to
# location given by $TMP

# The number of threads of RealignerTargetCreator can be controled with "-nt"
# SR0719r_404.MQ20.sorted.dedup.cv30.bam
#-fixMisencodedQuals \ #GATK commandline doc: Note that this argument should NEVER be used by default; you should only use it when you have confirmed that the quality scores in your data are not in the correct encoding.
java -Djava.io.tmpdir=$TMP -jar $GATKDIR/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-nt $CPU \
-R $REFFASTA \
-I $ACC.list \
-o AllSamples_forIndelRealigner.intervals

echo Realignment started
date

#-fixMisencodedQuals \
java -Djava.io.tmpdir=$TMP -jar $GATKDIR/GenomeAnalysisTK.jar \
-T IndelRealigner  \
-R $REFFASTA \
--targetIntervals AllSamples_forIndelRealigner.intervals \
-I $ACC.list \
--consensusDeterminationModel USE_READS  \
--nWayOut .realigned.bam

echo FixMateInformation started
date

# no list option for input, need loop
# Loop over each sample in bam.list, converting name to what is now the new ".realigned.bam" name

for BAMFILE in `cat bam.list`; do # egin loop by defining first SAMPLEFILE as the BAMFILE listed first in bam.list
SAMPLEBASE=`grep "$BAMFILE" bam.list | cut -d '.' -f1`    #find $BAMFILE in bam.list, then cut base name from first word before "."


java -Djava.io.tmpdir=$TMP -jar $PICARDDIR/picard.jar FixMateInformation \
INPUT=$SAMPLEBASE.MQ20.sorted.dedup.cv30.realigned.bam  \
OUTPUT=$SAMPLEBASE.MQ20.sorted.dedup.cv30.realigned.fixmate.bam \
SO=coordinate \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=$TMP

echo Run ended
date
done
```
---

start: Mon Feb 10 22:13:22 EST 2020

end: Wed Feb 12 20:58:20 EST 2020

Again, this wasn't very quick, for a total of 325 samples (I exclude 14 "rerun" samples in this step, with read counts ranging between 500K and 2M, and odd GC contents as thresholds. It took ~ 47 hours to finish. 
I also use a more strigent parameters to obtain a list of 311 samples only with reads > 2M

summary will be added here soon
