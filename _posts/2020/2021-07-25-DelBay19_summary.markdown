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

### Depth evaluation

Table 1。 Summary of the read depth distribution for each dataset. The last column is useful for Angsd -setMaxDepth setting.

| Data       | Mean | Deviation |  SD | Mean+3SD |
|------------|------|-----------|-----|----------|
| DelBay19   | 719  |   36664   | 191 |   1294   |
| DelBay20   | 190  |   2815    | 53  |   349    |
| DelBay19&20| 922  |   37462   | 194 |   1502   |
                 
DelBay19: dataset include DelBay19 challenge (n=97) and wild samples (n=234).             
DelBay20: dataset include DelBay20 challenge samples (n=101). 
DelBay19&20: dataset include DelBay19 and 20 samples (n=432)                  
       
### Sample coverage

[DelBay19_summary](https://github.com/hzz0024/HG_Code_Bay/blob/master/DelBay_summary/DelBay19_summary_final.xlsx)    

Mean depth (challenge): 0.78 (SD = 1.37)      

[DelBay20_summary](https://github.com/hzz0024/HG_Code_Bay/blob/master/DelBay_summary/DelBay20_summary_final.xlsx)  

Mean depth: 0.77 (SD = 1.34) 

### Global SNP calling using different datasets 

Table 5. Summary of SNPs used as a global SNP list (2032113 SNPs). Note that private and shared SNP number and the ratio are not the same between two datasets. I am working on other downsampling schemes (0.6 and 0.7x) to balance these numbers.

|                                          | Del19 (all samples)   | Del19 (only challenge)| Del20 (1x)              |
|------------------------------------------|-----------------------|-----------------------|-------------------------|
| Total number of sites analyzed           | 512960891             | 482669010             | 489141910               |
| Number of sites retained after filtering | 2322712               | 2140437               | 3600633                 |
| Private sites in each batch              | 290599 (12.51%)       | 379531 (17.73%)       | 1568520 (43.56%)        |
| Shared sites                             | 2032113 (87.49%)      | 1760906 (82.27%)      | 2032113 (56.44%)        |

### Relatedness

<img src="https://hzz0024.github.io/images/DelBay_adult/relatedness_wtoutlier.jpg" alt="img" width="800"/> 

### Diverstiy estimate

<img src="https://hzz0024.github.io/images/DelBay_adult/Diversity.jpg" alt="img" width="800"/>

### MDS

<img src="https://hzz0024.github.io/images/DelBay_adult/All_maf0.05_minq20_minmq25_pctind0.7_CV30_masked_noinvers_shared_sites.jpg" alt="img" width="800"/>

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

Here I am planning to update the reference database for SNP annotation and following enrichment analyses. A detailed protocal to build the whole reference dataset from scratch is [here](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=73#c), which includes three major steps: Diamond, InterproScan, and Blast2GO. However, to obtain this annotation database more efficiently (and more consistently with the current *Crassostrea virginica* genome manuscript), I just did some updates on the protein and gene description in original "Proestou_and_Sullivan_B2G_oyster_annotation" file. A full description about the this file is listed below,

- `  We used the NCBI protein sequence file, GCF00202275_2_c_virginica-3_0_protein.faa_gz, with the standard workflow and default parameters except we increased the GO weight parameter from 5 to 15.  We also used Kevin Johnson’s InterProScan results to supplement the blast results. The “Sequence Name” column is the protein product ID from the original NCBI protein sequence file used to run B2G.  The “Sequence Description” is from the B2G annotation run.  The columns that begin with “Annotation” contain the information on the GO terms that passed the annotation criteria AND were merged with Kevin’s InterProSCan inB2G (note: Part of the merging process is a validation step that removes redundant, more general functions based on the true path rule, only the most specific GO terms are assigned).  The columns that begin with “InterPro” contain only information imported from Kevin’s InterProScan.`

Below are some scripts for this update process:

1. python command to extrac the protein and gene description from faa and gtf files

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
# The results file is located in /Users/HG/Documents/HG/DelBay_adult/16_annotation/Data_base/deleted_GCF_002022765.2_C_virginica-3.0_feature_table.txt
```
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

Then combine everthing into a central database (done by Excel). There are some discrepancies between the updated and old (Dina's version) databases. For example, although the protein accession shared the same ID, the protein description is differed by **34747/60201(57.72%)**. I also observed **37 discrepancies** in gene ID (i.e. LOCXXXXXX). 

Currently I create a database with updated gene ID, protein and gene names (named GO_data_rm_dup.csv in /Users/HG/Documents/HG/DelBay_adult/16_annotation/Gene_annotation/). It includes 41928 records for SNP annotation.

2. python command to extract the gene description

A python script is created for annotation:

```python
import sys
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Input file.", default='XXX.bed')
parser.add_argument("-t", "--outlier", help="outlier file.", default='XXX.gene.txt')
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
        #lists.append((items[3], items[4], items[5], items[7].replace('\"',''), items[2])) # for protein annotation 
        lists.append((items[3], items[4], items[5], items[8].replace('\"',''), items[2])) # for gene annotation 

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

cat GEA_BF_20_Del19_FDR_2K.intersect # it only looks at the first three columns. Anything that overlaps with this regions will be recorded in the output file.
NC_035780.1 37977798    37979798    NC_035780.1_37978798    NC_035780.1 37977204    37979204    NC_035780.1_37978204    1406
NC_035783.1 14089954    14091954    NC_035783.1_14090954    NC_035783.1 14091394    14093394    NC_035783.1_14092394    560
NC_035784.1 41760936    41762936    NC_035784.1_41761936    NC_035784.1 41761900    41763900    NC_035784.1_41762900    1036
NC_035784.1 58317052    58319052    NC_035784.1_58318052    NC_035784.1 58316342    58318342    NC_035784.1_58317342    1290
NC_035786.1 8771371 8773371 NC_035786.1_8772371 NC_035786.1 87713718773371  NC_035786.1_8772371 2000
```

### Enrichment analysis using Gowinda

Gowinda is a multi-threaded Java application that allows unbiased analysis of gene set enrichment for Genome Wide Association Studies. Classical analysis of gene set (e.g.: Gene Ontology) enrichment assumes that all genes are sampled independently from each other with the same probability. These assumptions are violated in Genome Wide Association (GWA) studies since (i) longer genes typically have more SNPs resulting in a higher probability of being sampled and (ii) overlapping genes are sampled in clusters. Gowinda has been specifically designed to test for enrichment of gene sets in GWA studies. We show that Gene Ontology (GO) tests on GWA data could result in a substantial number of false positive GO terms. Permutation tests implemented in Gowinda eliminate these biases, but maintain sufficient power to detect enrichment of GO terms.

`Requirement`

Java 6 or higher.    
Furthermore the following input files are required    
a file containing the annotation of the genome in .gtf        
a gene set file, containing for every gene set (e.g.: GO category) a list of the associated gene IDs        
a file containing the total set of SNPs        
a file containing the candidate SNPs       


