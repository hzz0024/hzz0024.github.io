---
comments: true
title: DelBay19 poolseq analyses
date: '2020-08-11 12:00'
tags:
  - DelBay19
  - poolseq
  - ouliter
  - null model
  - WGS
categories:
  - WGS data analysis
---

### Reference about poolseq procedure

In [Evolutionary origins of genomic adaptations in an invasive copepod](https://www.nature.com/articles/s41559-020-1201-y), Stern and Lee (2020) used poolseq datasets to describe how biologial invaders, copepod, rapidly adapt to novel habitats. 

Their procedure for pool-seq data analysis:

1) raw read filtering and trimming -- BBDuk (BBTools package)    
2) genome mapping -- BEW-MEM v.0.7.17    
3) supplementary single-end reads mapping -- NextGenMap v.0.5.5    
4) duplicate reads removing -- Picard v. 2.18.27    
5) realign the regions around insertions and deletions -- GATK v3.8    
6) BAM file into mpileup format with Q>20 -- Samtools v.1.3.1    
7) sites within 3bp of an insertion or deletion were removed and the filtered mpileup conversion to sync format -- PoPoolation2    
8) bi-allelic SNP calling -- poolfstat v.1.0 (R), with parameters shown below,     
   MAF 0.05    
   at least four reads for a base call    
   a minimum of 20 and a maximum of 200 total read counts for all populations    

### Create a synchronized file

Given that I already have the generated bam files, I skipped the first six steps. I just merge all the bam file from one popualtion with samtools.

```sh
module load samtools/1.7
samtools merge HC.bam *.bam
samtools merge ARN.bam *.bam
samtools merge COH.bam *.bam
samtools merge SR.bam *.bam
samtools merge NB.bam *.bam
samtools merge CH.bam *.bam
samtools merge REF.bam *.bam
```

Note a warning message was shown in the log file: "PG tag "MarkDuplicates.9A" on read "NB551191:557:HVMVVBGXC:1:11112:14736:3819" encountered with no corresponding entry in header, tag lost. Unknown tags are only reported once per input file for each tag ID". The PG means the program information used for data producing. It seems there is some wrong with the header. For my current data analysis I will ignore this issue as it does not affect the output bam file. 

Synchronized files are the main input files for PoPoolation2. They basically contain the allele frequencies for every population at every base in the reference genome in a concise format. Note that the synchronized file format contains the allele frequencies after filtering for base quality. Here I created three population pair for trial analysis.

```sh 
#example in the tutorial file
samtools mpileup -B map/pop1.bam map/pop2.bam > p1_p2.mpileup
#for DelBay19 data
module load samtools/1.7
samtools mpileup -B CH.bam REF.bam > CH_REF.mpileup
samtools mpileup -B HC.bam NB.bam > HC_NB.mpileup
samtools mpileup -B COH.bam NB.bam > COH_NB.mpileup
# this step costs ~ 2 hours

# using mpileup2sync.pl to convert the mpileup into sync format
module load R/3.5.1
module load popoolation2/1201
perl /tools/popoolation2_1201/mpileup2sync.pl --fastq-type illumina --min-qual 20 --input CH_REF.mpileup --output CH_REF.sync

module load R/3.5.1
module load popoolation2/1201
perl /tools/popoolation2_1201/mpileup2sync.pl --fastq-type illumina --min-qual 20 --input COH_NB.mpileup --output COH_NB.sync


module load R/3.5.1
module load popoolation2/1201
perl /tools/popoolation2_1201/mpileup2sync.pl --fastq-type illumina --min-qual 20 --input HC_NB.mpileup --output HC_NB.sync

# for each of the formating run it costs ~ 8 hours.

# Note that Synchronizing the mpileup file is quite time consuming. To remove this bottleneck they implemented 'mpileup2sync' in Java multi-threading which is about 78x faster as the implementation in perl. However this option does not produce the right sync format file. May need further trial.

The output looks like this

NC_035780.1     25      N       0:0:7:0:0:0     0:0:1:0:0:0
NC_035780.1     26      N       12:0:0:0:0:0    2:0:0:0:0:0
NC_035780.1     27      N       0:13:0:0:0:0    0:2:0:0:0:0
NC_035780.1     28      N       14:0:0:0:0:0    2:0:0:0:0:0
NC_035780.1     29      N       0:0:12:0:0:0    0:0:2:0:0:0
NC_035780.1     30      N       1:0:0:10:0:0    0:0:0:2:0:0
NC_035780.1     31      N       0:14:0:0:0:0    0:1:0:0:0:0

col1: reference contig
col2: position within the refernce contig
col3: reference character
col4: allele frequencies of population number 1
col5: allele frequencies of population number 2
coln: allele frequencies of population number n

The allele frequencies are in the format A:T:C:G:N:del, i.e: count of bases 'A', count of bases 'T',... and deletion count in the end (character '*' in the mpileup)
```

### Calculate allele frequency differences

```sh
module load R/3.5.1
module load popoolation2/1201
perl /tools/popoolation2_1201/snp-frequency-diff.pl --input CH_REF.sync --output-prefix CH_REF1 --min-count 2 --min-coverage 5 --max-coverage 200
```

### Fisher's Exact Test for significance of absolute delta_p

The Fishers exact test can be used to test whether any differences in allele frequencies are statistically significant. At low coverages the absolute changes of allele frequencies or the Fst values may be strongly influenced by sampling effects, therefore the Fishers exact test may be used to identify significant changes in allele frequency.

```sh
module load R/3.5.1
module load popoolation2/1201
perl /tools/popoolation2_1201/fisher-test.pl --input CH_REF.sync --output CH_REF1.fet --min-count 2 --min-coverage 5 --max-coverage 200 --suppress-noninformative
```
