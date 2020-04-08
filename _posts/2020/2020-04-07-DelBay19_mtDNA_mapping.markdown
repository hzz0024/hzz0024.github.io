---
comments: true
title: DelBay19 mtDNA mapping
date: '2020-04-07 12:00'
tags:
  - samtools
  - Bowtie2
  - mapping
  - mitogenome
  - WGS
categories:
  - WGS data analysis
---

This post shows how to remove the reads that would map to the mitogenome. The detailed scripts & steps could be found in this [paper](https://research.bangor.ac.uk/portal/files/25243919/SNP_discovery_from_genomic_assemblies_Revised_191216.pdf), and its Github [site](https://github.com/PAMorin/SNPdiscovery). The major difference in my data running is that I did the mitogenome mapping by Bowtie2 but the paper used BWA for this purpose. 

---

- index the mtDNA genome

```sh
module load samtools/1.6 
module load bowtie2/2.3.5.1

# Index your reference using the faidx function in SAMtools  
samtools faidx cv30_mtDNA.fa

# Create index files needed for BWA read alignment to the reference   
bowtie2-build cv30_mtDNA.fa cv30_mtDNA >bt_mtDNA.log   
```
- start mapping & keeping unmapped reads

```sh
#!/bin/bash

# usage (from folder containing genome and trimmed folders): nohup ./Bowtie.sh REFERENCE REFBASENAME SAMPLELIST OUTPUTDIR >& mt_mapping.log

REFERENCE=$1 # path to reference fasta file
REFBASENAME=$2 # reference name to add to output files, e.g. cv30
SAMPLELIST=$3 # Filename and path to list of trimmed fastq files, e.g. samples.list
OUTPUTDIR=$4 # path to directory where output files are to be written

# Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST | cut -f4`; do #Begin loop by defining SAMPLEFILE as column 4 in $SAMPLELIST file, sample name + "_1_" is trimmed filename
SAMPLE=`grep "$SAMPLEFILE" $SAMPLELIST | cut -f4` #find $SAMPLEFILE in file $SAMPLELIST, then cut column 4 text >> SAMPLE (samplename)

#### MAPPING THE GENOMIC READS TO THE REFERENCE ####
# Extract relevant values from a table of seq ID, sample ID, library ID, and platform unit (here in columns 3, 4, and 2, respectively) for each sequenced library
SEQ_ID=`grep "${SAMPLEFILE}" $SAMPLELIST | cut -f1`
SAMPLE_ID=`grep "$SAMPLEFILE" $SAMPLELIST | cut -f4`
SAMPLE_SEQ_ID=$SAMPLE_ID'_'$SEQ_ID      # use when samples are sequenced in multiple lanes
LIB_ID=`grep "$SAMPLEFILE" $SAMPLELIST | cut -f2`
PU=`grep "$SAMPLEFILE" $SAMPLELIST | cut -f3`

# Map the paired-end reads: use --very-sensitive-local for reference transcriptome and --very-sensitive for reference genome
# --un-conc would output the unmapped reads in the current folder, with extension .1 (forward) and .2 (reverse)
bowtie2 -q --phred33 --very-sensitive --un-conc $SAMPLE -p 16 -I 0 -X 1500 --fr --rg-id $SAMPLE_SEQ_ID --rg SM:$SAMPLE_ID --rg LB:$LIB_ID --rg PU:$PU --rg PL:ILLUMINA -x /scratch/DelBay19/genome/cv30_mtDNA -1 $SAMPLE'_1_AdapterClipped_F_paired.fastq' -2 $SAMPLE'_1_AdapterClipped_R_paired.fastq' -S $SAMPLE'_Paired_bt2_'$REFBASENAME'.sam'

# ADD VIEW -Q 20 FILTER HERE
# Convert to bam file for storage
samtools view -q 20 -bS -F 4 $SAMPLE'_Paired_bt2_'$REFBASENAME'.sam' > $OUTPUTDIR'/'${SAMPLE}'_Paired_bt2_'${REFBASENAME}'.bam'
rm $SAMPLE'_Paired_bt2_'$REFBASENAME'.sam'

done
```
---

### Results

Overall, reads that could be mapped to the mitogenome contribute a very small portion of the trimmed data. The mitogenome mapping rates range from 0.00% to 3.78%, with an average value of 0.05%. The detailed summary data is [here](https://docs.google.com/spreadsheets/d/1VB4okrGUJDNhS5E9hsMUYVGbKtz8guPCs4SPiqKmqL4/edit?usp=sharing). 

Among them, two individuals showed abnormal high mapping rates, SR0719ch_002 (3,78%) and SR0719ch_327 (0.85%), which indicate potential contamination in these two samples. The low nuclear genome mapping rates of these two samples (SR0719ch_327 1.9% mapped and SR0719ch_002 SR002 60.89% mapped to the nuclear genome) relative to other samples (~80% mapping rate) supported this hypothesis.


The unmapped reads are used for the downstream nuclear DNA mapping and bam file processing. 





