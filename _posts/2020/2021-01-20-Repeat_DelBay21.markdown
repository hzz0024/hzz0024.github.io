---
comments: true
title: Repeat challenge sequencing analyses (DelBay20)
date: '2021-01-26 12:00'
tags:
  - DelBay
  - Challenge
  - WGS
  - QC
categories:
  - WGS data analysis
--- 

The dataset in this post comes from the repeated adult challenge experiment in summer, 2020. I am starting to trim the raw data and conduct read mapping using the same scripts. I also borrowed some of the code from Nina Lab [https://github.com/therkildsen-lab/genomic-data-analysis/blob/master/lcwgs_data_analysis.md](https://github.com/therkildsen-lab/genomic-data-analysis/blob/master/lcwgs_data_analysis.md)

Before data QC, we need to create two files to provide the detailed information of sequenced samples.

sample_list.txt: sample prefix part before *_1.fastq.gz* or *_2.fastq.gz*
sample_table.txt: a tab deliminated table with the following six columns, strictly in this order:

  - `prefix` the prefix of raw fastq file names

  - `lane_number` lane number; each sequencing lane or batch should be
    assigned a unique identifier

  - `seq_id` sequence ID; this variable is only relevant when different
    libraries were prepared out of the same sample and were run in the
    same lane (e.g. if you wanted to include a replicate). In this case,
    seq\_id should be used to distinguish these separate libraries. If
    you only have a single library prepared from each of your samples
    (even if you sequence that same library across multiple lanes), you
    can just put 1 for all samples in this column.

  - `sample_id` sample ID; a unique identifier for each individual
    sequenced

  - `population` population name; the population or other relevant
    grouping variable that the individual belongs to

  - `data_type` data type; there can only be two possible entries: `pe`
    (for paired-end data) or `se` (for single end data). We need this in
    the table because for some of our processing steps, the commands are
    slightly different for paired-end and single-end data.

These two files are located at            
[HG_Code_Bay/DelBay20/sample_list/fastq_list.txt](https://github.com/hzz0024/HG_Code_Bay/blob/master/DelBay20/sample_list/fastq_list.txt)                    
[HG_Code_Bay/DelBay20/sample_list/fastq_table.txt](https://github.com/hzz0024/HG_Code_Bay/blob/master/DelBay20/sample_list/fastq_table.txt)             

Besides, due to the long sample list, for each of the QC and trimming steps I created subsample lists (19) to reduce the running time. Therefore I have to write 19 scripts for each fo the step run. Thanks for for loop and echo functions, I made some script for that purpose. They are located at [HG_Code_Bay/DelBay20/script](https://github.com/hzz0024/HG_Code_Bay/tree/master/DelBay20/script), with wt_ as the prefix.

--- 

### Index-hopping issue

Index hopping occurs when dual indexes contain unmatched or unexpected pairs, and leading a small amount of contamination reads. This is a link to Youtube video explaining this phenomenon [https://www.youtube.com/watch?v=DR_8KbGGIhA&ab_channel=Illumina](https://www.youtube.com/watch?v=DR_8KbGGIhA&ab_channel=Illumina)

Here I checked the potential of index hopping in the sequencing data. 

```sh
zgrep "@" Cv5785_25_CKDL210000056-1a-AK11419-AK17135_HNMTKDSXY_L4_1.fq.gz | wc -l
17040229
zgrep "AATGGAGA+TGCACATA" Cv5785_25_CKDL210000056-1a-AK11419-AK17135_HNMTKDSXY_L4_1.fq.gz | wc -l
16855590
zgrep "@" Cv5785_25_CKDL210000056-1a-AK11419-AK17135_HNMTKDSXY_L4_1.fq.gz | grep -v "AATGGAGA+TGCACATA" | wc -l
184639
```

Code above calculate the number of reads with perfect match (16855590 of 17040229, 98.92%) after using unique dual indexes, and 184639 of 17040229 (1.08%) reads with 1 or 2 bp mismatch (corresponding the 2 bp error wiggle room used during demultiplexing). 

--- 

### Data process

1) First is to look at the fastqc reports:

```sh
#!/bin/sh
FASTQC=/programs/FastQC-0.11.8/fastqc
BASEDIR=/workdir/hz269/DelBay20
SAMPLELIST=$BASEDIR/sample_lists/fastq_list.txt
RAWFASTQSUFFIX1=_1.fq.gz # Suffix to raw fastq files. Use forward reads with paired-end data.

for SAMPLE in `cat $SAMPLELIST`; do

  $FASTQC $BASEDIR'/raw_fastq/'$SAMPLE$RAWFASTQSUFFIX1 -o $BASEDIR'/fastqc/'
done
```

2) Adapter clipping

Trimmomatic has lots of different filtering modules. Here I clip sequence that match to our adapter sequence (TruSeq3-PE-2.fa) and remove reads that end up being < 80bp after clipping.

```sh
cat 2_trim.sh
start=`date +%s`  ## date at start
BASEDIR=/workdir/hz269/DelBay20
TRIMMOMATIC=/programs/trimmomatic/trimmomatic-0.39.jar
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_test.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se. 
RAWFASTQDIR=$BASEDIR/raw_fastq/ # Path to raw fastq files. 
RAWFASTQSUFFIX1=_1.fq.gz # Suffix to raw fastq files. Use forward reads with paired-end data.
RAWFASTQSUFFIX2=_2.fq.gz # Suffix to raw fastq files. Use reverse reads with paired-end data. 
ADAPTERS=$BASEDIR/reference/NexteraPE-PE.fa # Path to a list of adapter/index sequences, copied from /programs/trimmomatic/adapters/

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do
    ## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
    SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
    POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
    SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
    LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
    SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely
    echo "Sample: $SAMPLE_UNIQ_ID"
    ## Extract data type from the sample table
    DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`
    
    ## The input and output path and file prefix
    RAWFASTQ_ID=$RAWFASTQDIR$SAMPLEFILE
    SAMPLEADAPT=$BASEDIR'/adapter_clipped/'$SAMPLE_UNIQ_ID
    
    ## Adapter clip the reads with Trimmomatic
    # The options for ILLUMINACLIP are: ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
    # The MINLENGTH drops the read if it is below the specified length in bp
    # For definitions of these options, see http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
    

    if [ $DATATYPE = pe ]; then
        java -jar $TRIMMOMATIC PE -threads 1 -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $RAWFASTQ_ID$RAWFASTQSUFFIX2 $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_f_unpaired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_unpaired.fastq.gz' 'ILLUMINACLIP:'$ADAPTERS':2:30:10:1:true MINLENGTH:80' 
    
    elif [ $DATATYPE = se ]; then
        java -jar $TRIMMOMATIC SE -threads 1 -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $SAMPLEADAPT'_adapter_clipped_se.fastq.gz' 'ILLUMINACLIP:'$ADAPTERS':2:30:10 MINLENGTH:40'
    fi
    
done

end=`date +%s`  ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 )) 
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"


    ## Adapter clip the reads with Trimmomatic
    # The options for ILLUMINACLIP are: ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>

# to run the script
for i in {1..19}; do
    nohup sh '2_trim_'$i'.sh' > '2_trim_'$i'.log' &
done
```

Runtime: 0:31:0 (hh:mm:ss) on sample Cv5785_25, with 34 M reads


3) Build reference index files

```sh
cat 3_build_ref.sh

BASEDIR=/workdir/hz269/DelBay20
PICARD=/programs/picard-tools-2.19.2/picard.jar
SAMTOOLS=/programs/samtools-1.11/bin/samtools
BOWTIEBUILD=/programs/bowtie2-2.3.5.1-linux-x86_64/bowtie2-build
REFERENCE=$BASEDIR/reference/CV30_masked.fasta   # This is a fasta file with the reference genome sequence we will map to
REFBASENAME="${REFERENCE%.*}"
$SAMTOOLS faidx $REFERENCE

java -jar $PICARD CreateSequenceDictionary R=$REFERENCE O=$REFBASENAME'.dict'
```

4) Trim polyg tails 

This is an option step in data process but could help trim the poly-G tails in the reads. PolyG is a common issue observed in Illumina NextSeq and NovaSeq series. It happens when a base with no light signal detected and called as "G", therefore poly-Gs is a sequencing artifacts.

```sh
#!/bin/bash
start=`date +%s`  ## date at start
## This script is used to quality filter and trim poly g tails. It can process both paired end and single end data.
BASEDIR=/workdir/hz269/DelBay19
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_1.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
FILTER=polyg # Type of filtering. Values can be: polyg (forced PolyG trimming only), quality (quality trimming, PolyG will be trimmed as well if processing NextSeq/NovaSeq data), or length (trim all reads to a maximum length)
THREAD=1 # Number of thread to use. Default is 10
FASTP=/programs/fastp-0.20.0/bin/fastp ## Path to the fastp program. The default path is /workdir/programs/fastp_0.19.7/fastp
MAXLENGTH=100 # Maximum length. This input is only relevant when FILTER=length, and its default value is 100.

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do

  ## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
  SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
  POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
  SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
  LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
  SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID
  echo "Sample: $SAMPLE_UNIQ_ID"
  ## Extract data type from the sample table
  DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`
  ## The input and output path and file prefix
  SAMPLEADAPT=$BASEDIR'/adapter_clipped/'$SAMPLE_UNIQ_ID
  SAMPLEQUAL=$BASEDIR'/qual_filtered/'$SAMPLE_UNIQ_ID

  ## Trim polyg tail or low quality tail with fastp.
  # --trim_poly_g forces polyg trimming, --cut_right enables cut_right quality trimming
  # -Q disables quality filter, -L disables length filter, -A disables adapter trimming
  # Go to https://github.com/OpenGene/fastp for more information
  if [ $DATATYPE = pe ]; then
    if [ $FILTER = polyg ]; then
      $FASTP --trim_poly_g --cut_right -L -A --thread $THREAD -i $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' -I $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_f_paired.fastq.gz' -O $SAMPLEQUAL'_adapter_clipped_qual_filtered_r_paired.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
    elif [ $FILTER = quality ]; then
      $FASTP -L -A --thread $THREAD -i $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' -I $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_f_paired.fastq.gz' -O $SAMPLEQUAL'_adapter_clipped_qual_filtered_r_paired.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
    elif [ $FILTER = length ]; then
      $FASTP --max_len1 $MAXLENGTH -Q -L -A --thread $THREAD -i $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' -I $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_f_paired.fastq.gz' -O $SAMPLEQUAL'_adapter_clipped_qual_filtered_r_paired.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
    fi
  fi
done

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
```

5) Map to the mtochondrial DNA

Map the paired-end reads to the mitochondrial DNA (ref_C_virginica-3.0_mtDNA.fasta.gz, here renamed as cv30_mtDNA) using Bowtie2. The remained reads will be mapped to nuclear DNA in step 6. I used --very-sensitive for reference genome mapping, and --un-conc OPTION to retain the PE reads that failed to align concordantly to the mtDNA genome

```sh
start=`date +%s` 
BASEDIR=/workdir/hz269/DelBay20
BOWTIE=/programs/bowtie2-2.3.5.1-linux-x86_64/bowtie2
SAMTOOLS=/programs/samtools-1.11/bin/samtools
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_test.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
FASTQDIR=$BASEDIR/qual_filtered/ # Path to the directory where fastq file are stored.
FASTQSUFFIX1=_adapter_clipped_qual_filtered_f_paired.fastq.gz  # Suffix to fastq files. Use forward reads with paired-end data.
FASTQSUFFIX2=_adapter_clipped_qual_filtered_r_paired.fastq.gz # Suffix to fastq files. Use reverse reads with paired-end data.
MAPPINGPRESET=very-sensitive # The pre-set option to use for mapping in bowtie2 (very-sensitive for end-to-end (global) mapping [typically used when we have a full genome reference], very-sensitive-local for partial read mapping that allows soft-clipping [typically used when mapping genomic reads to a transcriptome]
REFERENCE=$BASEDIR/reference/CV30_mtDNA.fasta # Path to reference fasta file and file name
REFNAME=CV30_mtDNA # Reference name to add to output files, e.g. gadMor2

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do

  ## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
  SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
  POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
  SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
  LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
  SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely

  ## Extract data type from the sample table
  DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`

  ## The input and output path and file prefix
  SAMPLETOMAP=$FASTQDIR$SAMPLE_UNIQ_ID
  SAMPLEBAM=$BASEDIR'/bam_mtDNA/'$SAMPLE_UNIQ_ID

  ## Define platform unit (PU), which is the lane number
  PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`

  ## Define reference base name
  REFBASENAME="${REFERENCE%.*}"

  ## Map reads to the reference
  echo $SAMPLE_UNIQ_ID

  # Map the mtDNA
  $BOWTIE -q --phred33 --$MAPPINGPRESET --un-conc-gz test -p 1 -I 0 -X 1500 --fr --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

done

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"

# to run the script
for i in {1..19}; do
    nohup sh '5_map_mtDNA_'$i'.sh' > '5_map_mtDNA_'$i'.log' &
done
```

6) Map to the reference, sort, and quality filter

```sh
start=`date +%s` 
BASEDIR=/workdir/hz269/DelBay20
BOWTIE=/programs/bowtie2-2.3.5.1-linux-x86_64/bowtie2
SAMTOOLS=/programs/samtools-1.11/bin/samtools
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_test.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
FASTQDIR=$BASEDIR/adapter_clipped/ # Path to the directory where fastq file are stored.
FASTQSUFFIX1=_adapter_clipped_f_paired.fastq.gz # Suffix to fastq files. Use forward reads with paired-end data.
FASTQSUFFIX2=_adapter_clipped_r_paired.fastq.gz # Suffix to fastq files. Use reverse reads with paired-end data.
MAPPINGPRESET=very-sensitive # The pre-set option to use for mapping in bowtie2 (very-sensitive for end-to-end (global) mapping [typically used when we have a full genome reference], very-sensitive-local for partial read mapping that allows soft-clipping [typically used when mapping genomic reads to a transcriptome]
REFERENCE=$BASEDIR/reference/CV30_masked.fasta # Path to reference fasta file and file name
REFNAME=CV30_masked # Reference name to add to output files, e.g. gadMor2

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do

  ## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
  SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
  POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
  SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
  LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
  SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely

  ## Extract data type from the sample table
  DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`

  ## The input and output path and file prefix
  SAMPLETOMAP=$FASTQDIR$SAMPLE_UNIQ_ID
  SAMPLEBAM=$BASEDIR'/bam/'$SAMPLE_UNIQ_ID

  ## Define platform unit (PU), which is the lane number
  PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`

  ## Define reference base name
  REFBASENAME="${REFERENCE%.*}"

  ## Map reads to the reference
  echo $SAMPLE_UNIQ_ID

  # Map the paired-end reads
  if [ $DATATYPE = pe ]; then
  # We ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (loose their mate so they become single-end), these can be mapped in a separate step and the resulting bam files merged with the paired-end mapped reads.
  $BOWTIE -q --phred33 --$MAPPINGPRESET -p 1 -I 0 -X 1500 --fr --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

  # Map the single-end reads
  elif [ $DATATYPE = se ]; then
  $BOWTIE -q --phred33 --$MAPPINGPRESET -p 1 --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -U $SAMPLETOMAP$FASTQSUFFIX1 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

  fi

  ## Convert to bam file for storage (including all the mapped reads)
  $SAMTOOLS view -bS -F 4 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam' > $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam'
  rm -f $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

  ## Filter the mapped reads (to onky retain reads with high mapping quality)
  # Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20)
  $SAMTOOLS view -h -q 20 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam' | $SAMTOOLS view -buS - | $SAMTOOLS sort -o $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'_minq20_sorted.bam'

done

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"

# to run the script
for i in {1..19}; do
    nohup sh '4_map_'$i'.sh' > '4_map_'$i'.log' &
done
```

Runtime: this process took so long. The rough eastimate for each individual is ~ 8 hours. Currently the sample list was divided into 19 part and ran in parallel, with ~ 10 samples per list. So the total time for the whole sequencing will take ~ 80 hours (3 days). 


7) Merge samples from different batches or lanes

8) Deduplicate and clip overlapping read pairs

```sh
#!/bin/bash
start=`date +%s` 
## This script is used to deduplicate bam files and clipped overlapping read pairs for paired end data. It can process both paired end and single end data.
BAMLIST=/workdir/hz269/DelBay_test/sample_lists/XXX # Path to a list of merged, deduplicated, and overlap clipped bam files. Full paths should be included. An example of such a bam list is /workdir/cod/greenland-cod/sample_lists/bam_list_1.tsv
SAMPLETABLE=$2 # Path to a sample table where the 1st column is the prefix of the MERGED bam files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The 5th column is population name and 6th column is the data type. An example of such a sample table is: /workdir/cod/greenland-cod/sample_lists/sample_table_merged.tsv
JAVA=java # Path to java
PICARD=/programs/picard-tools-2.19.2/picard.jar # Path to picard tools
BAMUTIL=/programs/bamUtil/bam # Path to bamUtil


## Loop over each sample
for SAMPLEBAM in `cat $BAMLIST`; do
  
  ## Extract the file name prefix for this sample
  SAMPLESEQID=`echo $SAMPLEBAM | sed 's/_bt2_.*//' | sed -e 's#.*/bam/\(\)#\1#'`
  SAMPLEPREFIX=`echo ${SAMPLEBAM%.bam}`

  ## Remove duplicates and print dupstat file
  # We used to be able to just specify picard.jar on the CBSU server, but now we need to specify the path and version
  $JAVA -Xmx60g -jar $PICARD MarkDuplicates I=$SAMPLEBAM O=$SAMPLEPREFIX'_dedup.bam' M=$SAMPLEPREFIX'_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
  
  ## Extract data type from the merged sample table
  DATATYPE=`grep -P "${SAMPLESEQID}\t" $SAMPLETABLE | cut -f 6`
  
  if [ $DATATYPE != se ]; then
    ## Clip overlapping paired end reads (only necessary for paired end data)
    $BAMUTIL clipOverlap --in $SAMPLEPREFIX'_dedup.bam' --out $SAMPLEPREFIX'_dedup_overlapclipped.bam' --stats
  fi
  
done

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
```

9) Indel realignment

```sh
#!/bin/bash
start=`date +%s`
## This script is used to quality filter and trim poly g tails. It can process both paired end and single end data. 
BAMLIST=BAMLIST=/workdir/hz269/DelBay_test/sample_lists/XXX # Path to a list of merged, deduplicated, and overlap clipped bam files. Full paths should be included. An example of such a bam list is /workdir/cod/greenland-cod/sample_lists/bam_list_1.tsv
BASEDIR=/workdir/hz269/DelBay_test # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories. An example for the Greenland cod data is: /workdir/cod/greenland-cod/
REFERENCE=$BASEDIR/reference/CV30_masked.fasta # Path to reference fasta file and file name, e.g /workdir/cod/reference_seqs/gadMor2.fasta
SAMTOOLS=/programs/samtools-1.11/bin/samtools # Path to samtools
GATK=${6:-/programs/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar} # Path to GATK

## Loop over each sample
for SAMPLEBAM in `cat $BAMLIST`; do

if [ -e $SAMPLEBAM'.bai' ]; then
  echo "the file already exists"
else
  ## Index bam files
  $SAMTOOLS index $SAMPLEBAM
fi

done

## Realign around in-dels
# This is done across all samples at once

## Use an older version of Java
export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH

## Create list of potential in-dels
if [ ! -f $BASEDIR'bam/all_samples_for_indel_realigner.intervals' ]; then
  java -Xmx40g -jar $GATK \
     -T RealignerTargetCreator \
     -R $REFERENCE \
     -I $BAMLIST \
     -o $BASEDIR'bam/all_samples_for_indel_realigner.intervals' \
     -drf BadMate
fi

## Run the indel realigner tool
java -Xmx40g -jar $GATK \
   -T IndelRealigner \
   -R $REFERENCE \
   -I $BAMLIST \
   -targetIntervals $BASEDIR'bam/all_samples_for_indel_realigner.intervals' \
   --consensusDeterminationModel USE_READS  \
   --nWayOut _realigned.bam

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
```

10) Estimate read depth

From Nina's github there are some code for depth estimation with samtools depth option. I previous used this method but confused about what coverage should I report. Perhaps I could adopt her method for sample coverage estimation. Link is [here](https://github.com/nt246/physalia-lcwgs/blob/main/day_1/markdowns/data_processing.md#estimate-read-depth-in-your-bam-files).

---

### Analyses summary

All processes were conducted on Hare CBSU server, with 48 cores, 256G RAM, and 144 TB disk space. This batch has 187 samples but divided into 19 groups for parallel running, with 10 samples in each group. To calculate total running time for each step just multiply 10 to average time (per sample). 

| Order | Step                          | Script Name    | Average time (per sample) | Note                                                 |
|-------|-------------------------------|----------------|---------------------------|------------------------------------------------------|
| 1     | FastQC                        | 1_fastqc.sh    |          25 mins          |                                                      |
| 2     | Adapter clipping              | 2_trim.sh      |          20 mins          |                                                      |
| 3     | Build reference               | 3_build_ref.sh |          < 5 mins         | Build both mtochondrial and nuclear reference genome |
| 4     | Trim polyg tails              | 4_polyg.sh     |          9 mins           |                                                      |
| 5     | Map to mtDNA reference        | 5_map_mtDNA.sh |          20 mins          |                                                      |
| 6     | Map to gDNA reference         | 6_map.sh       |          8 hours          |                                                      |
| 7     | Merge bam                     | 7_merge.sh     |                           | Only apply to DelBay19 data                          |
| 8     | Deduplicate and clip overlaps | 8_dup_clip.sh  |                           |                                                      |
| 9     | Indel realignment             | 9_realign.sh   |                           |                                                      |

---

### Example results for each step

1) adapter clipping: 

```sh
Sample: Cv104_CHR_1_4
Input Read Pairs: 22471817 Both Surviving: 22402666 (99.69%) Forward Only Surviving: 69058 (0.31%) Reverse Only Surviving: 91 (0.00%) Dro
pped: 2 (0.00%)
```

2) Trim polyg tails

```sh
Sample: Cv104_CHR_1_4

Read1 before filtering:
total reads: 22402666
total bases: 3360391893
Q20 bases: 3208650841(95.4844%)
Q30 bases: 2948308114(87.737%)

Read2 before filtering:
total reads: 22402666
total bases: 3360145453
Q20 bases: 3268116175(97.2612%)
Q30 bases: 3115230260(92.7112%)

Read1 after filtering:
total reads: 22301643
total bases: 3149733393
Q20 bases: 3050329153(96.844%)
Q30 bases: 2821412107(89.5762%)

Read2 aftering filtering:
total reads: 22301643
total bases: 3237927933
Q20 bases: 3191294433(98.5598%)
Q30 bases: 3057357416(94.4233%)

Filtering result:
reads passed filter: 44603286
reads failed due to low quality: 9376
reads failed due to too many N: 48

Duplication rate: 30.8626%
```

3) Map to mtDNA reference genome

```sh
Sample: Cv104_CHR_1_4

22301643 reads; of these:
  22301643 (100.00%) were paired; of these:
    22288067 (99.94%) aligned concordantly 0 times
    12027 (0.05%) aligned concordantly exactly 1 time
    1549 (0.01%) aligned concordantly >1 times
    ----
    22288067 pairs aligned concordantly 0 times; of these:
      142 (0.00%) aligned discordantly 1 time
    ----
    22287925 pairs aligned 0 times concordantly or discordantly; of these:
      44575850 mates make up the pairs; of these:
        44394533 (99.59%) aligned 0 times
        11127 (0.02%) aligned exactly 1 time
        170190 (0.38%) aligned >1 times
0.47% overall alignment rate
````
4) Map to masked reference genome

```sh
Sample: Cv104_CHR_1_4

22288067 reads; of these:
  22288067 (100.00%) were paired; of these:
    6919270 (31.04%) aligned concordantly 0 times
    9879276 (44.33%) aligned concordantly exactly 1 time
    5489521 (24.63%) aligned concordantly >1 times
    ----
    6919270 pairs aligned concordantly 0 times; of these:
      115700 (1.67%) aligned discordantly 1 time
    ----
    6803570 pairs aligned 0 times concordantly or discordantly; of these:
      13607140 mates make up the pairs; of these:
        9102058 (66.89%) aligned 0 times
        2360489 (17.35%) aligned exactly 1 time
        2144593 (15.76%) aligned >1 times
79.58% overall alignment rate
```

So far so good.