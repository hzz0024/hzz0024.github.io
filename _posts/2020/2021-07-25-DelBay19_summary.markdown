comments: true
title: DelBay 2019 summary
date: '2021-07-25 12:00'
tags:
  - DelBay19
  - PCA
  - Challenge
  - Wild 
  - WGS
  - QC
categories:
  - WGS data analysis
--- 

### Sample coverage

[DelBay19_summary](https://github.com/hzz0024/HG_Code_Bay/blob/master/DelBay_summary/DelBay19_summary_final.xlsx)    

| Data       | Mean coverage |   SD |Relealized coverage |  SD | 
|------------|---------------|------|--------------------|-----|
| DelBay19   | 0.91          | 0.38 |    2.49            | 0.68|

### Depth evaluation

Table 1。 Summary of the read depth distribution for each dataset. The last column is useful for Angsd -setMaxDepth setting.

| Data       | Mean | Deviation |  SD | Mean+3SD |
|------------|------|-----------|-----|----------|
| DelBay19   | 721  |   41946   | 205 |   1335   |
                 
DelBay19: dataset include DelBay19 challenge (n=97) and wild samples (n=234).

```sh
angsd command for global SNP calling
#!/bin/sh

###this script will work on bamfiles by population and calculate saf & maf
# maybe edit
target="Del19_final"
NB_CPU=24 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /local/workdir/hz269/DelBay19_angsd/01_scripts/01_config.sh
N_IND=$(wc -l $Del19_final | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "Ouput can be used for depth evaluation with all individuals listed in $Del19_final"
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is 70% of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

$angsd -P $NB_CPU \
 -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
 -doDepth 1 -doIBS 2 -makeMatrix 1 -doCov 1 $REGIONS \
 -maxDepth 2000 -dumpCounts 2 -anc $ANC -remove_bads 1 -rmTriallelic 1e-6 \
 -doPlink 2 -doGeno 4 -doPost 1 \
 -minMapQ 30 -minQ 20 -minInd $MIN_IND -minMaf $MIN_MAF \
 -setMinDepth 110 -setMaxDepth 1335 -SNP_pval 1e-6 -b $Del19_final \
 -out "/local/workdir/hz269/DelBay19_angsd/03_global/"$target"_maf"$MIN_MAF"_minq20_minmq30_pctind"$PERCENT_IND"_CV30_masked_noinvers"

 #main features
# -P nb of threads
# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2)
# -doMajorMinor 1 use the most frequent allele as major
# -anc provide a ancestral sequence = reference in our case
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

#main filters
#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ minimum quality of reads
#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%
#filter on allele frequency -minMaf, set to 0.05

#output

Total number of sites analyzed: 513108070
Number of sites retained after filtering: 2335739
[ALL done] cpu-time used =  221301.02 sec
[ALL done] walltime used =  53024.00 sec
```                            

### Relatedness

<img src="https://hzz0024.github.io/images/MDS/relatedness_wtoutlier.jpg" alt="img" width="800"/> 

### Diverstiy estimate

<img src="https://hzz0024.github.io/images/DelBay_adult/Diversity.jpg" alt="img" width="800"/>

### PCA & MDS

<img src="https://hzz0024.github.io/images/MDS/pca1.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/MDS/pca2.jpeg" alt="img" width="800"/>

### Combined Fisher's exact test

1) Update the depth for both Del19 (n=97) and Del20 (n=101) - done      
2) Rerun SNP calling for Del19 and Del20 - done    
3) Generate a global SNP list by extracting the common shared SNPs - done    
4) Calling allele frequency for each individual population using the global SNP list - done    
5) Perform Fisher’s exact tests - done
6) Examine common shared SNPs - done     
7) SNP annotation and functional gene exploration  

In the p-value combination step, because the weighted Z-test (also called ‘Stouffer’s method) is more power and more precision than does Fisher’s test ([Whitlock 2005](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2005.00917.x)), I used Z method to combine the p-values from the parallel tests. It favours symmetric rejection and is less sensitive to a single low p-value, requiring more consistently low p-values to yield a low combined p-value. 

The combinePValues function in the [scran R package](https://rdrr.io/bioc/scran/man/combinePValues.html) was used to perform Z method.

Results:

| Contrasts                  | No. outliers |
|----------------------------|--------------|
| REF19_CHR19_NB_HC          | 13           |
| REF19_CHR19_SR_HC          | 4            |
| Shared                     | 1            |
| REF19_SR_ARN_COH (control) | 0            |
| REF20_CHR20_NB_HC          | 7            |
| REF20_CHR20_SR_HC          | 0            |
| Shared                     | 0            |
| REF20_SR_ARN_COH (control) | 3            |
| REF19_REF20_CHR19_CHR20    | 696          |

The detailed SNP lists are [here](https://docs.google.com/spreadsheets/d/1hDH_lp_BQC9grGAK2vbZZBW-tzpiZCvMvajyMe-rYUo/edit?usp=sharing)

Now take a look at the allele frequency changes at these potential outliers.

- REF19_CHR19_NB_HC (13 outliers)

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_REF19_CHR19_NB_HC.jpg" alt="img" width="800"/>

One of them, SNP NC_035784.1_16552716 is shared between REF19_CHR19_NB_HC and REF19_CHR19_SR_HC. This SNPs is located in the gene Actin-depolymerizing factor 1-like (LOC111134891). Below is the delta_p patterns for this SNP

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_NC_035784.1_16552716.jpg" alt="img" width="800"/>

- REF20_CHR20_NB_HC (7 outliers)

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_REF20_CHR20_NB_HC.jpg" alt="img" width="800"/>

- REF19_REF20_CHR19_CHR20 (696 outlies)

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_control.jpg" alt="img" width="800"/>

Zoom in on some SNPs at chromosome 5

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_control_chr5.jpg" alt="img" width="800"/>

Zoom in on some SNPs at chromosome 1

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_control_chr1.jpg" alt="img" width="800"/>

### Single-generation selection (SGS) 

1) Determine the global theta value for each population (REF19, REF20, wild?) - REF19 done, the other two under running   
2) Determine -setMaxDepth for each population (CHR19, REF19, CHR20, REF20, HC, NB, SR) - done.      
3) Produce global SNP list for each of the population contrasts (CHR19 vs REF19, CHR20 vs REF20, HC vs. NB, HC vs SR) - done      
4) Generate allele frequency data for each of the population (with global SNP list, set as -site) - done      
5) Perform SGS tests on each of the contrast 

Table 6. Summary of SNP depth for Angsd -setMaxDepth settings. Note CHR19-REF19, CHR20-REF20, HC-SR, and HC-NB are contrasts used for SGS tests.

|     Contrasts       | Mean | Deviation |  SD | Mean+3SD |
|---------------------|------|-----------|-----|----------|
| CHR19-REF19         | 189  |   2966    | 54  |   352    |
| CHR20-REF20         | 433  |   33698   | 184 |   983    |
| HC-SR               | 216  |   4487    | 67  |   417    |
| HC_NB               | 210  |   3994    | 63  |   400    |

Note: p-value distribution is initially odd, probably due the usage of quantile-based p-value. After consulting with CSCU, now I switched to 2-sides p-value estimation with naive exterme delta-p counts.

p-value distribution for CHR19-REF19 contrast using 2-side p-value

<img src="https://hzz0024.github.io/images/DelBay_adult/ps_2side.jpeg" alt="img" width="800"/>

Table 6. Number of outliers identifed from SGS test for each population contrast (FDR < 0.05).

|     Contrasts       | Outliers |
|---------------------|----------|
| CHR19-REF19         | 3265     |
| CHR20-REF20         | 180      |              
| HC-SR               | 2444     |                     
| HC-NB               | 2559     |              

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_SGS.jpg" alt="img" width="800"/>

Shared SNPs:

CHR20-REF20 & HC_NB: NC_035782.1_62395550

CHR19-REF19 & HC_NB or HC_SR: 12 SNPS

CHR19-REF19 & HC_NB & HC_SR: NC_035786.1_8772371

### Probabilistic Random Forest

Following the paper by Reis et al. 2018. [Probabilistic Random Forest: A machine learning algorithm for noisy datasets](https://arxiv.org/pdf/1811.05994.pdf). I am trying to incorporate the genotype likelihood into random forest test.

The initial trial is performed on 3006 outliers SNPs identified from SGS CHR19-REF19 contrasts (total 2032113 SNPS).

<img src="https://hzz0024.github.io/images/DelBay_adult/PRF_accuracy.jpg" alt="img" width="800"/>

- shared SNPs among repeat runs (1000, 500, 100, 50 are top SNPs based on importance)

1000
NC_035780.1_50001008 NC_035782.1_45223277 NC_035784.1_10248759 *NC_035784.1_16552962* NC_035784.1_18676441 NC_035784.1_83138159 NC_035786.1_33698717

500
NC_035780.1_50001008 NC_035782.1_45223277 *NC_035784.1_16552962* NC_035784.1_18676441 NC_035786.1_33698717

100
*NC_035784.1_16552962*

50
*NC_035784.1_16552962*

- Zoom-in for NC_035784.1_16552962 (red) and NC_035784.1_16552716 (lightgreen), which both located in Actin-depolymerizing factor 1-like (LOC111134891).

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_PRF.jpg" alt="img" width="800"/>

### Genotype-environment association

- RDA

```sh
# prepare the test data
head -n 1001 by_pop_0.05_pctind0.7_maxdepth3.mafs.rda.bak > by_pop_0.05_pctind0.7_maxdepth3.mafs.rda # only take the first 1000 SNPs

````

### SNP Annotation

NCBI recently released some tools to retrieve gene sequence and metadata. Based on NCBI, currently ***Crassostrea virginica*** genome has 39,493 genes and pseudogenes, which includes 34,596 protein-coding genes, 4,230 non-coding genes, and 667 pseudogenes. All these gene_id can be obtained from the [GCF_002022765.2_C_virginica-3.0_feature_table.txt.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_feature_table.txt.gz) file. 

There are some discrepancies between the updated and old (Dina's version) databases. For example, although the protein accession shared the same ID, the protein description is differed by **34747/60201(57.72%)**. I also observed **37 discrepancies** in gene_id (note the gene_id means the numbers after LOC). 

Currently I create a database with updated gene ID and gene names (named Gene_annotation_all.csv in /Users/HG/Documents/HG/DelBay_adult/16_annotation/Gene_annotation//Users/HG/Documents/HG/DelBay_adult/16_annotation/Gene_annotation_all.csv). It includes 39,493 records.

Below are the details for Gene_annotation_all.csv creation:

Download the NCBI tools

NCBI Datasets command line tools are datasets and dataformat [link](https://www.ncbi.nlm.nih.gov/datasets/docs/command-line/).

- datasets to download biological sequence data across all domains of life from NCBI.

- dataformat to convert metadata from JSON Lines format to other formats.

```bash

1. R script to create the bash command
setwd("/Users/HG/Documents/HG/DelBay_adult/16_annotation/Data_base_format/test")
file = 'dataset_39493.txt' 
d = read.delim(file, header = T, sep='')
dt = as.vector(t(d))
# split into 800 entries per run
max <- 800
x <- seq_along(dt)
dd <- split(dt, ceiling(x/max))
dd
length(dd[[1]])
paste(as.character(dd[[1]]), collapse=",")

vv = c()
for(i in 1:50) {
  v = paste(as.character(dd[[i]]), collapse=",")
  v_c = paste0("./datasets download gene gene-id ",v, " --filename ", i, ".zip" )
  vv = c(vv, v_c)
}

write.table(vv, "./process.sh", row.names=F, quote=F, sep="\n")

# format
kk = c()
for(i in 1:50) {
  k = paste0("./dataformat tsv gene --package ",i, ".zip --fields chromosomes,common-name,description,ensembl-geneids,gene-id,gene-type,genomic-range-range-orientation,genomic-range-range-start,genomic-range-range-stop,genomic-region-gene-range-range-start,genomic-region-gene-range-range-stop,genomic-region-genomic-region-type,tax-id > ", i, ".tsv" )
  kk = c(kk, k)
}
write.table(kk, "./format.sh", row.names=F, quote=F, sep="\n")

2. ./process.sh and ./format.sh

3. cat *.tsv > total.tsv | mv total.tsv Gene_annotation_all.csv

```

2. python command to extract the gene description

A python script is created for annotation:

```python
import sys
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Input file.", default='GO_data_format.csv')
parser.add_argument("-t", "--outlier", help="outlier file.", default='CS_HC_noinvers.sliding.zfst.outlier.merged.igv')
parser.add_argument("-o", "--output", help="Output", default='GO_target.csv')
args = parser.parse_args()

infile = args.input
outlier_file = args.outlier
outfile = args.output

import csv
from io import StringIO
# read outliers
lists = []
with open(infile, 'r') as f:
    HEADER = False
    for l in f:
        if not HEADER:
            HEADER = True
            continue
        reader = csv.reader(StringIO(l), delimiter=',')
        for r in reader:
            items = r
        #items = l.split(',') 
        lists.append((items[1], items[2], items[3], items[4].replace('\"',''), items[5])) # for gene annotation 

with open(outlier_file, 'r') as f, open(outfile, 'w') as w:

    for l in f:
        OVERLAP = False
        #if not HEADER:
        #    HEADER = True
        #    continue
        items = l.split()
        chr = items[0]
        st = int(items[1])
        ed = int(items[2])
        annos = []
        geneIDs = []
        for _chr, _st, _ed, _anno, _GeneID in lists:
            _st = int(_st)
            _ed = int(_ed)
            if _chr == chr and st <= _ed and ed >= _st:
                OVERLAP = True
                annos.append(_anno)
                geneIDs.append(_GeneID)

        if OVERLAP:
             annos = list(set(annos))
             geneIDs = list(set(geneIDs))
             outline = chr + '\t' + str(st) + '\t' + str(ed) + '\t' + str(len(geneIDs)) + '\t' + ';'.join(annos) + '\t' + ';'.join(geneIDs)
             w.write(outline + '\n') 
        else:
             outline = chr + '\t' + str(st) + '\t' + str(ed) + '\t' + 'No mapped genes'
             w.write(outline + '\n') 


usage example
python3 extract_gene_V2.py -i GO_data_rm_dup.csv -t GEA_BF_20_Del19_FDR_2K.intersect -o GEA_BF_20_Del19_FDR_2K.intersect.gene.txt        

# input example
cat GEA_BF_20_Del19_FDR_2K.intersect # it only looks at the first three columns. Anything that overlaps with this regions will be recorded in the output file.
NC_035780.1 37977798    37979798    NC_035780.1_37978798    NC_035780.1 37977204    37979204    NC_035780.1_37978204    1406
NC_035783.1 14089954    14091954    NC_035783.1_14090954    NC_035783.1 14091394    14093394    NC_035783.1_14092394    560
NC_035784.1 41760936    41762936    NC_035784.1_41761936    NC_035784.1 41761900    41763900    NC_035784.1_41762900    1036
NC_035784.1 58317052    58319052    NC_035784.1_58318052    NC_035784.1 58316342    58318342    NC_035784.1_58317342    1290
NC_035786.1 8771371 8773371 NC_035786.1_8772371 NC_035786.1 87713718773371  NC_035786.1_8772371 2000

# output example
cat GEA_BF_20_Del19_FDR_2K.intersect.gene.txt
NC_035780.1 37977798    37979798    1   disks large-associated protein 4-like   111123126
NC_035783.1 14089954    14091954    1   tubulin polyglutamylase ttll6-like  111130353
NC_035784.1 41760936    41762936    2   carbohydrate sulfotransferase 3-like;brevican core protein-like 111135454;111135457
NC_035784.1 58317052    58319052    1   cholecystokinin receptor-like   111132944
NC_035786.1 8771371 8773371 1   neuron navigator 2-like 111103516
```

### Enrichment analysis using Gowinda

#### 1. Data description

A detailed protocal to build the whole GO reference dataset from scratch is [here](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=73#c), which includes three major steps: Diamond, InterproScan, and Blast2GO. 

Here I only did some updates on the protein and gene description in original "Proestou_and_Sullivan_B2G_oyster_annotation" file. This GO reference is build by Dina's team and a full description about the this file is listed below,


- `  We used the NCBI protein sequence file, GCF00202275_2_c_virginica-3_0_protein.faa_gz, with the standard workflow and default parameters except we increased the GO weight parameter from 5 to 15.  We also used Kevin Johnson’s InterProScan results to supplement the blast results. The “Sequence Name” column is the protein product ID from the original NCBI protein sequence file used to run B2G.  The “Sequence Description” is from the B2G annotation run.  The columns that begin with “Annotation” contain the information on the GO terms that passed the annotation criteria AND were merged with Kevin’s InterProSCan inB2G (note: Part of the merging process is a validation step that removes redundant, more general functions based on the true path rule, only the most specific GO terms are assigned).  The columns that begin with “InterPro” contain only information imported from Kevin’s InterProScan.`

#### 2. Gowinda

Gowinda is a multi-threaded Java application that allows unbiased analysis of gene set enrichment for Genome Wide Association Studies. Classical analysis of gene set (e.g.: Gene Ontology) enrichment assumes that all genes are sampled independently from each other with the same probability. These assumptions are violated in Genome Wide Association (GWA) studies since (i) longer genes typically have more SNPs resulting in a higher probability of being sampled and (ii) overlapping genes are sampled in clusters. Gowinda has been specifically designed to test for enrichment of gene sets in GWA studies. We show that Gene Ontology (GO) tests on GWA data could result in a substantial number of false positive GO terms. Permutation tests implemented in Gowinda eliminate these biases, but maintain sufficient power to detect enrichment of GO terms.

`Requirement`

Java 6 or higher.    
Furthermore the following input files are required    

- a file containing the annotation of the genome in .gtf   

```sh
# process the gtf file
# gene id
cat GCF_002022765.2_C_virginica-3.0_genomic.gtf |awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | tail -n +5 > tmp1.txt
cat GCF_002022765.2_C_virginica-3.0_genomic.gtf |awk -F "\t" '{print $9}'|awk -F " " '{print $1" "$2}' | tail -n +5 > tmp2.txt
paste -d "\t" tmp1.txt tmp2.txt > tmp.gtf
grep -v "NC_007175.2" tmp.gtf > annotation.gtf
rm tmp*
# delete any residual content in the annotation.gtf file (e.g. the ### at the bottom)
```

- a gene set file, containing for every gene set (e.g.: GO category) a list of the associated gene IDs  

```python
filename = 'GO_all_need_format.txt'
outfile = 'GO_format.txt'

GO_map = {}
GO_term_map = {}
i = 0
HEADER = True
with open(filename, 'r') as f:
    for l in f:
        if HEADER:
            HEADER = False
            continue
        i += 1
        items = l.strip().split('\t')
        if len(items) < 3:
            continue
        assert len(items) == 3
        Gene_ID = items[0]
        GO_IDs = items[1].strip('\"').split(';')
        GO_terms = items[2].strip('\"').split(';')
        assert len(GO_IDs) == len(GO_terms)
        for GO_ID, GO_term in zip(GO_IDs, GO_terms):
            if GO_ID in GO_term_map:
                assert GO_term == GO_term_map[GO_ID]
            else:
                GO_term_map[GO_ID] = GO_term

        if GO_ID in GO_map:
            GO_map[GO_ID].append(Gene_ID)
        else:
            GO_map[GO_ID] = [Gene_ID]

GO_IDs = list(GO_map.keys())
GO_IDs = sorted(GO_IDs, key=lambda x: int(x.split(':')[-1]))
with open(outfile, 'w') as w:
    for GO_ID in GO_IDs:
        w.write(GO_ID + '\t' + GO_term_map[GO_ID] + '\t' + ' '.join(GO_map[GO_ID]) + '\n')   

head GO_all_need_format.txt

Gene_ID GO_ID   GO_term
LOC111138521    GO:0004731;GO:0009116   purine-nucleoside phosphorylase activity;nucleoside metabolic process
LOC111099031
LOC111099031
LOC111099031
LOC111099032    GO:0004930;GO:0007186;GO:0016021    G protein-coupled receptor activity;G protein-coupled receptor signaling pathway;integral component of membrane
LOC111099034
LOC111099035
LOC111099035
LOC111099033    GO:0016021;GO:0022857;GO:0055085    integral component of membrane;transmembrane transporter activity;transmembrane transport

head GO_format.txt
GO:0000062  fatty-acyl-CoA binding  LOC111130021 LOC111130071
GO:0000123  histone acetyltransferase complex   LOC111127833
GO:0000166  nucleotide binding  LOC111105732 LOC111105732 LOC111105954 LOC111106157 LOC111106157 LOC111107229 LOC111108558 LOC111108558 LOC111109807 LOC111111648 LOC111113560 LOC111116119 LOC111116119 LOC111116494 LOC111118596 LOC111118742 LOC111120214 LOC111120886 LOC111123184 LOC111123885 LOC111123885 LOC111124351 LOC111125312 LOC111125389 LOC111125389 LOC111127527 LOC111129267 LOC111129267 LOC111129337 LOC111129337 LOC111132062 LOC111132328 LOC111132330 LOC111132455 LOC111132974 LOC111133887 LOC111133906 LOC111134584 LOC111135002 LOC111136494 LOC111136669 LOC111136669 LOC111136669 LOC111136669 LOC111136669 LOC111137206 LOC111137464 LOC111137734
GO:0000184  nuclear-transcribed mRNA catabolic process, nonsense-mediated decay LOC111136887 LOC111136887 LOC111138487
GO:0000226  microtubule cytoskeleton organization   LOC111132121 LOC111132123
GO:0000290  deadenylation-dependent decapping of nuclear-transcribed mRNA   LOC111100321 LOC111100321 LOC111101096 LOC111101096
GO:0000350  generation of catalytic spliceosome for second transesterification step LOC111121655
GO:0000398  mRNA splicing, via spliceosome  LOC111131404
GO:0000422  autophagy of mitochondrion  LOC111118945 LOC111118945
GO:0000462  maturation of SSU-rRNA from tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)    LOC111110205
```

- a file containing the total set of SNPs

```sh
# process the total SNPs
total_file=CHR19_all_minq20_minmq30_CV30_masked.mafs
cat $total_file |awk -F "\t" '{print $1"\t"$2}' | tail -n +2 > total_snps.txt
```

- a file containing the candidate SNPs       

```sh
# Process the candidate SNPs
cand_file=salinity2_env2.txt_candidate_SNP_2_sd.txt
cat $cand_file |awk '{print $2}'|awk -F "_" '{print $1"_"$2"\t"$3}' | tail -n +2 > cand_snps.txt
```

Note: need to replace the chromosome name with numbers

```sh
for i in annotation.gtf cand_snps.txt total_snps.txt; do
sed -i .bak 's/NC_035780.1/1/g;s/NC_035781.1/2/g;s/NC_035782.1/3/g;s/NC_035783.1/4/g;s/NC_035784.1/5/g;s/NC_035785.1/6/g;s/NC_035786.1/7/g;s/NC_035787.1/8/g;s/NC_035788.1/9/g;s/NC_035789.1/10/g' $i
done

./sed_chr.sh
rm *.bak
```

The files are located in /Users/HG/Documents/HG/DelBay_adult/16_annotation/GO_enrichment. Now run the script below

```sh
# Basic analysis
java -Xmx4g -jar ./Gowinda-1.12.jar --snp-file total_snps.txt --candidate-snp-file cand_snps.txt --gene-set-file GO_format.txt --annotation-file annotation.gtf --simulations 100000 --min-significance 1 --gene-definition gene --threads 8 --output-file results_gene_gene.txt --mode gene --min-genes 1
# Including regulatory regions
java -Xmx4g -jar ./Gowinda-1.12.jar --snp-file total_snps.txt --candidate-snp-file cand_snps.txt --gene-set-file GO_format.txt --annotation-file annotation.gtf --simulations 100000 --min-significance 1 --gene-definition updownstream5000 --threads 8 --output-file results_gene_5000ud.txt --mode gene --min-genes 1
# high resolution GO term enrichment
java -Xmx4g -jar ./Gowinda-1.12.jar --snp-file total_snps.txt --candidate-snp-file cand_snps.txt --gene-set-file GO_format.txt --annotation-file annotation.gtf --simulations 100000 --min-significance 1 --gene-definition updownstream5000 --threads 8 --output-file results_snp_5000ud.txt --mode snp --min-genes 1
```

Check [https://code.google.com/archive/p/gowinda/wikis/Manual.wiki](https://code.google.com/archive/p/gowinda/wikis/Manual.wiki) for data preparation and result interpretation [https://sourceforge.net/p/gowinda/wiki/Tutorial/](https://sourceforge.net/p/gowinda/wiki/Tutorial/)

#### 3. Enrichment tests

- outliers in RDA analysis

RDA has identifed a total of 

---

Other useful command:

1. python command to extract the protein description from C_virginica-3.0_feature_table file. 

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_feature_table.txt.gz
#subset the feature table (done by Excel), the file name is /Users/HG/Documents/HG/DelBay_adult/16_annotation/Data_base/GCF_002022765.2_C_virginica-3.0_feature_table.txt

# below is python3 command
# to delete the isoform X from orthologs
filename = 'GCF_002022765.2_C_virginica-3.0_feature_table_edit.txt'
outfile = 'isoform_deleted_GCF_002022765.2_C_virginica-3.0_feature_table.txt'

HEADER = True
with open(filename, 'r') as f, open(outfile, 'w') as w:
    for l in f:
        if HEADER:
            HEADER = False
            continue
        items = l.strip().split('\t') 
        protein_name = items[5]
        if 'isoform' in protein_name:
            items[5] = protein_name.split('isoform')[0].strip()
            w.write('\t'.join(items) + '\n')
        else:
            w.write(l)
# The results file is located in /Users/HG/Documents/HG/DelBay_adult/16_annotation/Data_base/isoform_deleted_GCF_002022765.2_C_virginica-3.0_feature_table.txt
```

1. python command to extract the protein description from faa files. 

```python
# download the protein fasta: 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_protein.faa.gz
# below is python3 command
filename = 'GCF_002022765.2_C_virginica-3.0_protein.faa' # this is the most recent protein fasta
outfile = 'GCF_002022765.2_C_virginica-3.0_protein_annotation.txt'
with open(filename, 'r') as f, open(outfile, 'w') as w:
    i = 0
    for l in f:
        if l.startswith('>'):
            i += 1
            l = l.strip().strip('>')
            items = l.split('[')[0].split(' ', 1)
            w.write('\t'.join(items) + '\n')
        else:
            continue
# note there is still some error in the processed file, contents from GCF_002022765.2_C_virginica-3.0_feature_table.txt.gz is more accurate.
```
