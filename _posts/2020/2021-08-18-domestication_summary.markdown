---
comments: true
title: 600K Domestication Project summary
date: '2021-08-18 12:00'
tags:
  - 600K
  - SNP array
  - oyster
  - domestication 
  - WGS
  - summary
categories:
  - WGS data analysis
--- 

### 1. VCF process

Filtering parameters: --maf 0.05, --max-missing 0.7, --chr 1-10

Number of SNPs in the original vcf: 300,446  
Number of total samples: 842    
Number of SNPs after filtering: 159,849      
Number of total samples: 514

#### 1.1 Vcf file evaluation

Figure 1. SNP summary for vcf file. From top to bottom: Proportion of missing data per SNP; Proportion of missing data per individual; distribution of minor allele frequencies

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/vcf_summary.jpg" alt="img" width="800"/>

Our vcf data has a very promising profile for downstream analyses - 1) clearly most individuals have a call at almost every site; 2) the proportion of missing data per individual is low; 3) it is clear that a large number of variants have low frequency alleles.

### 2. Fst and ZFst



### 3. Diverstiy estimate



### 4. PCA & MDS

<img src="https://hzz0024.github.io/images/MDS/pca1.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/MDS/pca2.jpeg" alt="img" width="800"/>

### 5. Combined Fisher's exact test

### 6. Random Forest

Following the paper by Reis et al. 2018. [Probabilistic Random Forest: A machine learning algorithm for noisy datasets](https://arxiv.org/pdf/1811.05994.pdf). I am trying to incorporate the genotype likelihood into random forest test.

#### delta_p patterns

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

First, let's test whether the candidate SNPs show any enrichment for any GO category when counting every gene only once (--mode gene). Furthermore we want a SNP to be associated with a gene if the SNP is overlapping with exons or introns (--gene-definition gene). To speed up the tutorial let's use 8 threads for the simulations and perform 100.000 simulations.

```
Output headers:

column 1: the GO term
column 2: on the average this number of genes are found per simulation for the given GO category. In --mode gene every gene is only counted once whereas in --mode snp a single gene may be counted several times dependent on the SNP
column 3: using the candidate SNPs this number of genes was found for the given GO category. In --mode gene every gene is only counted once whereas in --mode snp a single gene may be counted several times dependent on the SNP
column 4: p-value (uncorrected for multiple testing)
column 5: FDR (p-value after adjustment for multiple testing)
column 6: the number of genes (uniq) found for the given GO category
column 7: the number of genes that could at most be found for the given GO category, i.e.: genes of the given GO category that have an corresponding entry in the annotation file and contain at least one SNP
column 8: total number of genes for the given GO category in the GO association file
column 9: description of the given GO term
column 10: comma separated list of the gene_ids found for the given GO category
```

```sh
# Basic analysis 
java -Xmx4g -jar ./Gowinda-1.12.jar --snp-file total_snps.txt --candidate-snp-file cand_snps.txt --gene-set-file GO_format.txt --annotation-file annotation.gtf --simulations 100000 --min-significance 1 --gene-definition gene --threads 8 --output-file results_gene_gene.txt --mode gene --min-genes 1

# below is cut with p-value < 0.1
GO:0006265  0.325   3   0.0009800000    0.1825800000    3   3   3   DNA topological change  loc111135726,loc111113980,loc111125129
GO:1902751  0.007   1   0.0067000000    0.7606791176    1   1   7   positive regulation of cell cycle G2/M phase transition loc111132091
GO:1903087  0.010   1   0.0103200000    0.7606791176    1   1   2   mitotic spindle pole body duplication   loc111128092
GO:0016746  1.070   4   0.0160600000    0.7606791176    4   18  32  transferase activity, transferring acyl groups  loc111130055,loc111117093,loc111103965,loc111137103
GO:0005509  7.822   14  0.0175400000    0.7606791176    14  108 154 calcium ion binding loc111133254,loc111126968,loc111135892,loc111135463,loc111132196,loc111125455,loc111101038,loc111129884,loc111129988,loc111119733,loc111119729,loc111129799,loc111128063,loc111135845
GO:0045947  0.018   1   0.0182500000    0.7606791176    1   1   1   negative regulation of translational initiation loc111121086
GO:0097546  0.022   1   0.0223300000    0.7606791176    1   1   1   ciliary base    loc111120088
GO:0061542  0.024   1   0.0243100000    0.7606791176    1   1   2   3-demethylubiquinol-n 3-O-methyltransferase activity    loc111104443
GO:0106105  0.025   1   0.0249200000    0.7606791176    1   1   1   Ala-tRNA(Thr) hydrolase activity    loc111129835
GO:0048382  0.026   1   0.0255300000    0.7606791176    1   1   2   mesendoderm development loc111102517
GO:0044057  0.026   1   0.0260200000    0.7606791176    1   1   1   regulation of system process    loc111133412
GO:0032039  0.294   2   0.0269800000    0.7606791176    2   3   5   integrator complex  loc111123220,loc111133970
GO:0097433  0.029   1   0.0287000000    0.7606791176    1   1   1   dense body  loc111135917
GO:0004197  0.344   2   0.0331700000    0.7606791176    2   3   3   cysteine-type endopeptidase activity    loc111124470,loc111117344
GO:0017148  0.034   1   0.0344300000    0.7606791176    1   1   1   negative regulation of translation  loc111121303
GO:0050265  0.035   1   0.0347800000    0.7606791176    1   1   1   RNA uridylyltransferase activity    loc111099079
GO:0005686  0.036   1   0.0355000000    0.7606791176    1   1   1   U2 snRNP    loc111136243
GO:0047631  0.036   1   0.0355300000    0.7606791176    1   1   2   ADP-ribose diphosphatase activity   loc111136071
GO:0030127  0.337   2   0.0365200000    0.7606791176    2   4   4   COPII vesicle coat  loc111118731,loc111120467
GO:0007275  0.908   3   0.0367300000    0.7606791176    3   6   7   multicellular organism development  loc111133485,loc111114450,loc111118493
GO:0006091  0.038   1   0.0372500000    0.7606791176    1   2   2   generation of precursor metabolites and energy  loc111125475
GO:0007165  7.984   13  0.0415200000    0.7606791176    13  80  151 signal transduction loc111124558,loc111114535,loc111136106,loc111115210,loc111103265,loc111119258,loc111121590,loc111124631,loc111137522,loc111134242,loc111125571,loc111118443,loc111130167
GO:0017056  0.374   2   0.0444200000    0.7606791176    2   4   5   structural constituent of nuclear pore  loc111135936,loc111125592
GO:0050936  0.046   1   0.0463900000    0.7606791176    1   1   1   xanthophore differentiation loc111127266
GO:0010468  0.410   2   0.0464300000    0.7606791176    2   3   3   regulation of gene expression   loc111120529,loc111119382
GO:1990090  0.047   1   0.0470000000    0.7606791176    1   1   2   cellular response to nerve growth factor stimulus   loc111135421
GO:0071310  0.047   1   0.0474400000    0.7606791176    1   1   1   cellular response to organic substance  loc111136514
GO:0043142  0.048   1   0.0478200000    0.7606791176    1   1   1   single-stranded DNA helicase activity   loc111132235
GO:0006383  0.048   1   0.0479900000    0.7606791176    1   1   2   transcription by RNA polymerase III loc111130603
GO:0070161  0.048   1   0.0481900000    0.7606791176    1   1   2   anchoring junction  loc111110132
GO:0070300  0.049   1   0.0493700000    0.7606791176    1   1   1   phosphatidic acid binding   loc111136102
GO:0006782  0.052   1   0.0510900000    0.7606791176    1   2   2   protoporphyrinogen IX biosynthetic process  loc111118050
GO:0004364  0.053   1   0.0521100000    0.7606791176    1   4   7   glutathione transferase activity    loc111106298
GO:2000543  0.052   1   0.0521400000    0.7606791176    1   1   2   positive regulation of gastrulation loc111121370
GO:0032456  0.057   1   0.0567900000    0.8272765714    1   1   1   endocytic recycling loc111131029
GO:0008146  0.932   3   0.0604000000    0.8529856757    3   18  26  sulfotransferase activity   loc111110499,loc111123814,loc111130175
GO:0016573  0.061   1   0.0607200000    0.8529856757    1   3   15  histone acetylation loc111123211
GO:0047429  0.065   1   0.0644500000    0.9006436842    1   2   3   nucleoside-triphosphate diphosphatase activity  loc111118714
GO:0042256  0.069   1   0.0687200000    0.9068089362    1   1   1   mature ribosome assembly    loc111138435
GO:0046854  0.466   2   0.0694800000    0.9068089362    2   7   11  phosphatidylinositol phosphorylation    loc111112578,loc111128712
GO:0006807  0.071   1   0.0707300000    0.9068089362    1   1   1   nitrogen compound metabolic process loc111108142
GO:0071986  0.072   1   0.0715700000    0.9068089362    1   2   6   Ragulator complex   loc111119487
GO:0046662  0.073   1   0.0726700000    0.9068089362    1   1   1   regulation of oviposition   loc111136938
GO:0043086  0.075   1   0.0749600000    0.9068089362    1   1   5   negative regulation of catalytic activity   loc111111017
GO:0061077  0.077   1   0.0757000000    0.9068089362    1   3   4   chaperone-mediated protein folding  loc111127664
GO:0090148  0.076   1   0.0759600000    0.9068089362    1   1   1   membrane fission    loc111130294
GO:1903311  0.076   1   0.0764300000    0.9068089362    1   1   1   regulation of mRNA metabolic process    loc111131549
GO:0070266  0.078   1   0.0781200000    0.9137477083    1   1   1   necroptotic process loc111130239
GO:0060255  0.084   1   0.0824900000    0.9415548077    1   2   2   regulation of macromolecule metabolic process   loc111103921
GO:0098845  0.084   1   0.0836900000    0.9415548077    1   1   1   postsynaptic endosome   loc111108538
GO:0097431  0.085   1   0.0851800000    0.9415548077    1   1   1   mitotic spindle pole    loc111121417
GO:0030246  3.859   7   0.0851900000    0.9415548077    7   77  131 carbohydrate binding    loc111129223,loc111130470,loc111109234,loc111119285,loc111116829,loc111129951,loc111129106
GO:0070449  0.089   1   0.0872000000    0.9463856757    1   2   3   elongin complex loc111125674
GO:0046474  0.089   1   0.0891300000    0.9463856757    1   1   1   glycerophospholipid biosynthetic process    loc111099557
GO:0007264  1.857   4   0.0904200000    0.9463856757    4   13  21  small GTPase mediated signal transduction   loc111125960,loc111124657,loc111108426,loc111130373
GO:0032501  0.094   1   0.0920500000    0.9463856757    1   2   4   multicellular organismal process    loc111109628
```

In the previous analysis only SNPs overlapping with exons or introns were associated with a gene. In order to account for regulatory regions, let's also include 5000 bp upstream and 5000 bp downstream of a gene (--gene-definition updownstream5000), for the next analysis.

```sh
# Including regulatory regions up and downstream 5000 bp
java -Xmx4g -jar ./Gowinda-1.12.jar --snp-file total_snps.txt --candidate-snp-file cand_snps.txt --gene-set-file GO_format.txt --annotation-file annotation.gtf --simulations 100000 --min-significance 1 --gene-definition updownstream5000 --threads 8 --output-file results_gene_5000ud.txt --mode gene --min-genes 1

GO:0006265  0.459   3   0.0032300000    0.4396125000    3   3   3   DNA topological change  loc111135726,loc111113980,loc111125129
GO:0005509  12.472  22  0.0055800000    0.4396125000    22  113 154 calcium ion binding loc111133254,loc111126968,loc111100577,loc111135892,loc111135463,loc111132196,loc111125455,loc111101038,loc111129884,loc111121949,loc111112623,loc111134696,loc111130740,loc111138432,loc111129988,loc111119733,loc111119729,loc111129799,loc111129797,loc111128063,loc111120363,loc111135845
GO:0070449  0.163   2   0.0057600000    0.4396125000    2   2   3   elongin complex loc111125674,loc111132184
GO:0016746  1.950   6   0.0069300000    0.4396125000    6   18  32  transferase activity, transferring acyl groups  loc111130055,loc111117093,loc111103965,loc111100498,loc111101391,loc111137103
GO:0061077  0.145   2   0.0102400000    0.5702320000    2   4   4   chaperone-mediated protein folding  loc111127664,loc111127012
GO:0071819  0.271   2   0.0177100000    0.8024875000    2   2   2   DUBm complex    loc111103564,loc111133642
GO:0042030  0.019   1   0.0190700000    0.8024875000    1   1   1   ATPase inhibitor activity   loc111116974
GO:0071949  1.648   5   0.0202200000    0.8024875000    5   24  34  FAD binding loc111130282,loc111113982,loc111104633,loc111136227,loc111131016
GO:0030170  3.656   8   0.0257800000    0.8024875000    8   36  50  pyridoxal phosphate binding loc111130191,loc111111022,loc111138313,loc111128100,loc111119686,loc111137345,loc111130121,loc111128099
GO:0046854  0.798   3   0.0329700000    0.8024875000    3   8   11  phosphatidylinositol phosphorylation    loc111112578,loc111128712,loc111118897
GO:1902751  0.034   1   0.0336000000    0.8024875000    1   1   7   positive regulation of cell cycle G2/M phase transition loc111132091
GO:0048382  0.034   1   0.0340100000    0.8024875000    1   1   2   mesendoderm development loc111102517
GO:0008180  0.329   2   0.0345700000    0.8024875000    2   4   5   COP9 signalosome    loc111136748,loc111127772
GO:0050265  0.037   1   0.0368000000    0.8024875000    1   1   1   RNA uridylyltransferase activity    loc111099079
GO:0005730  0.809   3   0.0376400000    0.8024875000    3   9   13  nucleolus   loc111120892,loc111129049,loc111138077
GO:0106105  0.039   1   0.0393400000    0.8024875000    1   1   1   Ala-tRNA(Thr) hydrolase activity    loc111129835
GO:0016592  2.069   5   0.0419500000    0.8024875000    5   16  23  mediator complex    loc111125864,loc111137835,loc111133969,loc111137901,loc111131411
GO:0098556  0.048   1   0.0479100000    0.8024875000    1   1   1   cytoplasmic side of rough endoplasmic reticulum membrane    loc111125240
GO:0034450  0.399   2   0.0481300000    0.8024875000    2   3   4   ubiquitin-ubiquitin ligase activity loc111135880,loc111121932
GO:0042602  0.050   1   0.0495000000    0.8024875000    1   1   1   riboflavin reductase (NADPH) activity   loc111121274
GO:1903087  0.050   1   0.0497300000    0.8024875000    1   1   2   mitotic spindle pole body duplication   loc111128092
GO:0008146  1.408   4   0.0517000000    0.8024875000    4   19  26  sulfotransferase activity   loc111114753,loc111110499,loc111123814,loc111130175
GO:0043142  0.053   1   0.0533800000    0.8024875000    1   1   1   single-stranded DNA helicase activity   loc111132235
GO:0098519  0.054   1   0.0537500000    0.8024875000    1   1   1   nucleotide phosphatase activity, acting on free nucleotides loc111128130
GO:0008305  0.805   3   0.0541200000    0.8024875000    3   9   9   integrin complex    loc111113860,loc111109550,loc111112736
GO:0045947  0.055   1   0.0552300000    0.8024875000    1   1   1   negative regulation of translational initiation loc111121086
GO:0017148  0.055   1   0.0553600000    0.8024875000    1   1   1   negative regulation of translation  loc111121303
GO:0060840  0.055   1   0.0554000000    0.8024875000    1   1   3   artery development  loc111113798
GO:0070300  0.056   1   0.0562900000    0.8024875000    1   1   1   phosphatidic acid binding   loc111136102
GO:2000772  0.058   1   0.0578200000    0.8024875000    1   1   1   regulation of cellular senescence   loc111133475
GO:0031204  0.058   1   0.0583700000    0.8024875000    1   1   2   posttranslational protein targeting to membrane, translocation  loc111101095
GO:0032039  0.453   2   0.0583800000    0.8024875000    2   3   5   integrator complex  loc111123220,loc111133970
GO:0006098  0.060   1   0.0597600000    0.8024875000    1   1   2   pentose-phosphate shunt loc111105262
GO:0032367  0.383   2   0.0639000000    0.8024875000    2   4   4   intracellular cholesterol transport loc111129681,loc111130084
GO:0006996  0.065   1   0.0654700000    0.8024875000    1   1   1   organelle organization  loc111135980
GO:0043461  0.066   1   0.0659300000    0.8024875000    1   1   1   proton-transporting ATP synthase complex assembly   loc111121651
GO:0007165  11.689  17  0.0660500000    0.8024875000    17  95  151 signal transduction loc111124558,loc111114535,loc111125225,loc111136106,loc111115210,loc111103265,loc111119258,loc111121590,loc111124631,loc111137522,loc111134242,loc111125571,loc111103359,loc111103360,loc111138472,loc111118443,loc111130167
GO:2000543  0.068   1   0.0680700000    0.8024875000    1   1   2   positive regulation of gastrulation loc111121370
GO:0030127  0.475   2   0.0686500000    0.8024875000    2   4   4   COPII vesicle coat  loc111118731,loc111120467
GO:0048034  0.069   1   0.0693900000    0.8024875000    1   1   1   heme O biosynthetic process loc111117954
GO:0044459  0.072   1   0.0708800000    0.8024875000    1   2   2   obsolete plasma membrane part   loc111116151
GO:1904888  0.478   2   0.0715900000    0.8024875000    2   4   5   cranial skeletal system development loc111107844,loc111123841
GO:0005686  0.073   1   0.0728300000    0.8024875000    1   1   1   U2 snRNP    loc111136243
GO:0030042  0.077   1   0.0751200000    0.8024875000    1   2   2   actin filament depolymerization loc111128666
GO:0032434  0.076   1   0.0761300000    0.8024875000    1   1   1   regulation of proteasomal ubiquitin-dependent protein catabolic process loc111128109
GO:0010468  0.522   2   0.0767900000    0.8024875000    2   3   3   regulation of gene expression   loc111120529,loc111119382
GO:0032549  0.080   1   0.0785100000    0.8024875000    1   3   4   ribonucleoside binding  loc111112947
GO:0004197  0.525   2   0.0786500000    0.8024875000    2   3   3   cysteine-type endopeptidase activity    loc111124470,loc111117344
GO:0097546  0.079   1   0.0787600000    0.8024875000    1   1   1   ciliary base    loc111120088
GO:0050545  0.080   1   0.0798000000    0.8024875000    1   1   1   sulfopyruvate decarboxylase activity    loc111114715
GO:0010466  0.468   2   0.0802600000    0.8024875000    2   6   13  negative regulation of peptidase activity   loc111117429,loc111116198
GO:0002682  0.081   1   0.0806700000    0.8024875000    1   1   1   regulation of immune system process loc111120716
GO:0007264  2.474   5   0.0811400000    0.8024875000    5   15  21  small GTPase mediated signal transduction   loc111125960,loc111109273,loc111124657,loc111108426,loc111130373
GO:0047631  0.081   1   0.0813500000    0.8024875000    1   1   2   ADP-ribose diphosphatase activity   loc111136071
GO:0046662  0.083   1   0.0826200000    0.8024875000    1   1   1   regulation of oviposition   loc111136938
GO:0006782  0.085   1   0.0834700000    0.8024875000    1   2   2   protoporphyrinogen IX biosynthetic process  loc111118050
GO:0070181  0.084   1   0.0836700000    0.8024875000    1   1   2   small ribosomal subunit rRNA binding    loc111120557
GO:0017056  0.521   2   0.0836700000    0.8024875000    2   4   5   structural constituent of nuclear pore  loc111135936,loc111125592
GO:0051539  1.739   4   0.0847000000    0.8024875000    4   17  29  4 iron, 4 sulfur cluster binding    loc111121293,loc111117101,loc111118241,loc111135616
GO:0007275  1.198   3   0.0872900000    0.8024875000    3   6   7   multicellular organism development  loc111133485,loc111114450,loc111118493
GO:0006995  0.088   1   0.0878800000    0.8024875000    1   1   1   cellular response to nitrogen starvation    loc111134873
GO:0016757  6.247   10  0.0884200000    0.8024875000    10  61  103 transferase activity, transferring glycosyl groups  loc111110121,loc111133477,loc111134324,loc111125433,loc111128714,loc111125440,loc111137074,loc111099105,loc111120548,loc111134973
GO:0050688  0.089   1   0.0886100000    0.8024875000    1   2   3   regulation of defense response to virus loc111129814
GO:0030942  0.091   1   0.0886800000    0.8024875000    1   2   2   endoplasmic reticulum signal peptide binding    loc111130615
GO:0050936  0.090   1   0.0899200000    0.8207281538    1   1   1   xanthophore differentiation loc111127266
GO:0071310  0.091   1   0.0913900000    0.8254635714    1   1   1   cellular response to organic substance  loc111136514
GO:0035097  0.607   2   0.0923600000    0.8254635714    2   2   2   histone methyltransferase complex   loc111099693,loc111125290
GO:0044057  0.095   1   0.0952400000    0.8254635714    1   1   1   regulation of system process    loc111133412
GO:0046777  0.573   2   0.0958700000    0.8254635714    2   4   4   protein autophosphorylation loc111124479,loc111127183
GO:1990259  0.096   1   0.0961500000    0.8254635714    1   1   2   histone-glutamine methyltransferase activity    loc111122686
GO:0032021  0.098   1   0.0979900000    0.8254635714    1   1   1   NELF complex    loc111125986
GO:0035725  2.615   5   0.0980100000    0.8254635714    5   19  29  sodium ion transmembrane transport  loc111119600,loc111129263,loc111128402,loc111120535,loc111111628
GO:2000601  0.100   1   0.0996100000    0.8254635714    1   1   2   positive regulation of Arp2/3 complex-mediated actin nucleation loc111108074
```

As we are getting quite desperate and we know that linkage disequilibrium is decaying quite rapidly in our organism we want to perform a higher resolution analysis in which genes may be counted several times, dependent on the number of candidate SNPs (--mode snp). For example a gene containing 5 candidate SNPs is counted 5 times (with the simulations as well as with the candidate SNPs).

```sh
# high resolution GO term enrichment
java -Xmx4g -jar ./Gowinda-1.12.jar --snp-file total_snps.txt --candidate-snp-file cand_snps.txt --gene-set-file GO_format.txt --annotation-file annotation.gtf --simulations 100000 --min-significance 1 --gene-definition updownstream5000 --threads 8 --output-file results_snp_5000ud.txt --mode snp --min-genes 1

GO:0051601  0.191   3   0.0010700000    0.3618766667    1   1   1   exocyst localization    loc111127094
GO:0042602  0.052   2   0.0013800000    0.3618766667    1   1   1   riboflavin reductase (NADPH) activity   loc111121274
GO:0016746  2.201   8   0.0017700000    0.3618766667    6   18  32  transferase activity, transferring acyl groups  loc111130055,loc111117093,loc111103965,loc111100498,loc111101391,loc111137103
GO:0043890  0.351   3   0.0052400000    0.8079807895    1   1   2   N-acetylgalactosamine-6-sulfatase activity  loc111134371
GO:0097433  0.117   2   0.0062900000    0.8079807895    1   1   1   dense body  loc111135917
GO:0060828  0.456   3   0.0111100000    0.8079807895    1   2   3   regulation of canonical Wnt signaling pathway   loc111126518
GO:0070266  0.161   2   0.0111200000    0.8079807895    1   1   1   necroptotic process loc111130239
GO:1901505  0.166   2   0.0120800000    0.8079807895    1   1   2   carbohydrate derivative transmembrane transporter activity  loc111131860
GO:0070449  0.174   2   0.0134000000    0.8079807895    2   2   3   elongin complex loc111125674,loc111132184
GO:0031224  0.179   2   0.0145300000    0.8079807895    1   1   1   intrinsic component of membrane loc111133518
GO:0061077  0.155   2   0.0148500000    0.8079807895    2   4   4   chaperone-mediated protein folding  loc111127664,loc111127012
GO:0006265  0.511   3   0.0153300000    0.8079807895    3   3   3   DNA topological change  loc111135726,loc111113980,loc111125129
GO:0042030  0.019   1   0.0192200000    0.8079807895    1   1   1   ATPase inhibitor activity   loc111116974
GO:0017056  0.575   3   0.0203000000    0.8079807895    2   4   5   structural constituent of nuclear pore  loc111135936,loc111125592
GO:1990380  0.219   2   0.0204600000    0.8079807895    1   1   1   Lys48-specific deubiquitinase activity  loc111129798
GO:0043236  0.218   2   0.0207100000    0.8079807895    1   1   1   laminin binding loc111110868
GO:0010468  0.592   3   0.0222900000    0.8079807895    2   3   3   regulation of gene expression   loc111120529,loc111119382
GO:2000146  0.242   2   0.0238100000    0.8079807895    1   1   1   negative regulation of cell motility    loc111121005
GO:0000422  0.240   2   0.0245200000    0.8079807895    1   1   1   autophagy of mitochondrion  loc111118945
GO:0030170  4.064   9   0.0264200000    0.8079807895    8   36  50  pyridoxal phosphate binding loc111130191,loc111111022,loc111138313,loc111128100,loc111119686,loc111137345,loc111130121,loc111128099
GO:0005509  14.424  23  0.0279200000    0.8079807895    22  113 154 calcium ion binding loc111133254,loc111126968,loc111100577,loc111135892,loc111135463,loc111132196,loc111125455,loc111101038,loc111129884,loc111121949,loc111112623,loc111134696,loc111130740,loc111138432,loc111129988,loc111119733,loc111119729,loc111129799,loc111129797,loc111128063,loc111120363,loc111135845
GO:0046777  0.654   3   0.0285200000    0.8079807895    2   4   4   protein autophosphorylation loc111124479,loc111127183
GO:0044238  0.681   3   0.0311300000    0.8079807895    2   9   24  primary metabolic process   loc111119606,loc111120493
GO:0006471  0.683   3   0.0329200000    0.8079807895    2   6   11  protein ADP-ribosylation    loc111100311,loc111135118
GO:0048382  0.034   1   0.0337100000    0.8079807895    1   1   2   mesendoderm development loc111102517
GO:1902751  0.035   1   0.0346600000    0.8079807895    1   1   7   positive regulation of cell cycle G2/M phase transition loc111132091
GO:0071819  0.298   2   0.0361400000    0.8079807895    2   2   2   DUBm complex    loc111103564,loc111133642
GO:0050265  0.039   1   0.0377700000    0.8079807895    1   1   1   RNA uridylyltransferase activity    loc111099079
GO:0071949  1.798   5   0.0381400000    0.8079807895    5   24  34  FAD binding loc111130282,loc111113982,loc111104633,loc111136227,loc111131016
GO:0006364  0.740   3   0.0390600000    0.8079807895    2   4   5   rRNA processing loc111128153,loc111130503
GO:0048266  0.311   2   0.0395500000    0.8079807895    1   1   1   behavioral response to pain loc111128509
GO:0106105  0.041   1   0.0402500000    0.8079807895    1   1   1   Ala-tRNA(Thr) hydrolase activity    loc111129835
GO:0051377  0.317   2   0.0413600000    0.8079807895    1   3   3   mannose-ethanolamine phosphotransferase activity    loc111119347
GO:0046907  0.764   3   0.0419400000    0.8079807895    2   7   12  intracellular transport loc111135919,loc111129530
GO:1901068  0.329   2   0.0437900000    0.8079807895    1   2   3   guanosine-containing compound metabolic process loc111125192
GO:0006464  2.085   6   0.0444500000    0.8079807895    4   9   10  cellular protein modification process   loc111135074,loc111135085,loc111125930,loc111130353
GO:0072546  0.336   2   0.0448900000    0.8079807895    1   5   7   ER membrane protein complex loc111136484
GO:0008295  0.337   2   0.0456300000    0.8079807895    1   2   2   spermidine biosynthetic process loc111103445
GO:0098556  0.049   1   0.0474600000    0.8117573077    1   1   1   cytoplasmic side of rough endoplasmic reticulum membrane    loc111125240
GO:0008180  0.358   2   0.0501400000    0.8117573077    2   4   5   COP9 signalosome    loc111136748,loc111127772
GO:1903087  0.052   1   0.0515500000    0.8117573077    1   1   2   mitotic spindle pole body duplication   loc111128092
GO:0043142  0.055   1   0.0539600000    0.8117573077    1   1   1   single-stranded DNA helicase activity   loc111132235
GO:0060840  0.058   1   0.0566000000    0.8117573077    1   1   3   artery development  loc111113798
GO:0098519  0.059   1   0.0572000000    0.8117573077    1   1   1   nucleotide phosphatase activity, acting on free nucleotides loc111128130
GO:0045947  0.059   1   0.0576100000    0.8117573077    1   1   1   negative regulation of translational initiation loc111121086
GO:0017148  0.060   1   0.0580100000    0.8117573077    1   1   1   negative regulation of translation  loc111121303
GO:2000772  0.060   1   0.0580600000    0.8117573077    1   1   1   regulation of cellular senescence   loc111133475
GO:0070300  0.060   1   0.0586600000    0.8117573077    1   1   1   phosphatidic acid binding   loc111136102
GO:0031204  0.061   1   0.0587300000    0.8117573077    1   1   2   posttranslational protein targeting to membrane, translocation  loc111101095
GO:0006098  0.062   1   0.0603800000    0.8117573077    1   1   2   pentose-phosphate shunt loc111105262
GO:0005730  0.887   3   0.0605000000    0.8117573077    3   9   13  nucleolus   loc111120892,loc111129049,loc111138077
GO:0046854  0.889   3   0.0608500000    0.8117573077    3   8   11  phosphatidylinositol phosphorylation    loc111112578,loc111128712,loc111118897
GO:0006996  0.068   1   0.0654200000    0.8377959184    1   1   1   organelle organization  loc111135980
GO:0043066  0.420   2   0.0670500000    0.8377959184    1   5   9   negative regulation of apoptotic process    loc111107995
GO:0043461  0.072   1   0.0690800000    0.8377959184    1   1   1   proton-transporting ATP synthase complex assembly   loc111121651
GO:2000543  0.072   1   0.0696700000    0.8377959184    1   1   2   positive regulation of gastrulation loc111121370
GO:0007165  13.721  20  0.0711300000    0.8377959184    17  95  151 signal transduction loc111124558,loc111114535,loc111125225,loc111136106,loc111115210,loc111103265,loc111119258,loc111121590,loc111124631,loc111137522,loc111134242,loc111125571,loc111103359,loc111103360,loc111138472,loc111118443,loc111130167
GO:0048034  0.074   1   0.0711800000    0.8377959184    1   1   1   heme O biosynthetic process loc111117954
GO:0034450  0.438   2   0.0722800000    0.8377959184    2   3   4   ubiquitin-ubiquitin ligase activity loc111135880,loc111121932
GO:0044459  0.075   1   0.0728500000    0.8377959184    1   2   2   obsolete plasma membrane part   loc111116151
GO:0005686  0.077   1   0.0737400000    0.8377959184    1   1   1   U2 snRNP    loc111136243
GO:0008146  1.526   4   0.0745800000    0.8377959184    4   19  26  sulfotransferase activity   loc111114753,loc111110499,loc111123814,loc111130175
GO:0030042  0.081   1   0.0774900000    0.8377959184    1   2   2   actin filament depolymerization loc111128666
GO:0008305  0.893   3   0.0787600000    0.8377959184    3   9   9   integrin complex    loc111113860,loc111109550,loc111112736
GO:0032434  0.083   1   0.0794300000    0.8377959184    1   1   1   regulation of proteasomal ubiquitin-dependent protein catabolic process loc111128109
GO:0097546  0.083   1   0.0798400000    0.8377959184    1   1   1   ciliary base    loc111120088
GO:0050545  0.085   1   0.0818600000    0.8377959184    1   1   1   sulfopyruvate decarboxylase activity    loc111114715
GO:0032549  0.086   1   0.0822600000    0.8377959184    1   3   4   ribonucleoside binding  loc111112947
GO:0002682  0.087   1   0.0831700000    0.8377959184    1   1   1   regulation of immune system process loc111120716
GO:0046983  9.135   14  0.0835300000    0.8377959184    12  75  108 protein dimerization activity   loc111111083,loc111107401,loc111105492,loc111105253,loc111119882,loc111133713,loc111131330,loc111129050,loc111102494,loc111119796,loc111131265,loc111134789
GO:0047631  0.088   1   0.0843300000    0.8377959184    1   1   2   ADP-ribose diphosphatase activity   loc111136071
GO:0032367  0.425   2   0.0858000000    0.8377959184    2   4   4   intracellular cholesterol transport loc111129681,loc111130084
GO:0046662  0.090   1   0.0861500000    0.8377959184    1   1   1   regulation of oviposition   loc111136938
GO:0007264  2.969   6   0.0864500000    0.8377959184    5   15  21  small GTPase mediated signal transduction   loc111125960,loc111109273,loc111124657,loc111108426,loc111130373
GO:0016592  2.333   5   0.0873300000    0.8377959184    5   16  23  mediator complex    loc111125864,loc111137835,loc111133969,loc111137901,loc111131411
GO:0006782  0.091   1   0.0874400000    0.8377959184    1   2   2   protoporphyrinogen IX biosynthetic process  loc111118050
GO:0070181  0.093   1   0.0886000000    0.8377959184    1   1   2   small ribosomal subunit rRNA binding    loc111120557
GO:0030942  0.093   1   0.0889900000    0.8377959184    1   2   2   endoplasmic reticulum signal peptide binding    loc111130615
GO:0050688  0.093   1   0.0892800000    0.8377959184    1   2   3   regulation of defense response to virus loc111129814
GO:0006995  0.094   1   0.0898700000    0.8377959184    1   1   1   cellular response to nitrogen starvation    loc111134873
GO:0050936  0.097   1   0.0923700000    0.8377959184    1   1   1   xanthophore differentiation loc111127266
GO:0071310  0.098   1   0.0932100000    0.8377959184    1   1   1   cellular response to organic substance  loc111136514
GO:0032039  0.508   2   0.0933900000    0.8377959184    2   3   5   integrator complex  loc111123220,loc111133970
GO:0008289  1.670   4   0.0943800000    0.8377959184    1   15  21  lipid binding   loc111125856
GO:0016757  6.915   11  0.0951600000    0.8377959184    10  61  103 transferase activity, transferring glycosyl groups  loc111110121,loc111133477,loc111134324,loc111125433,loc111128714,loc111125440,loc111137074,loc111099105,loc111120548,loc111134973
GO:1904888  0.522   2   0.0979200000    0.8377959184    2   4   5   cranial skeletal system development loc111107844,loc111123841
GO:0030127  0.525   2   0.0984600000    0.8377959184    2   4   4   COPII vesicle coat  loc111118731,loc111120467
GO:0016887  4.621   8   0.0984700000    0.8377959184    5   28  37  ATPase activity loc111122859,loc111119223,loc111123982,loc111103560,loc111118893
GO:0044057  0.104   1   0.0991500000    0.8377959184    1   1   1   regulation of system process    loc111133412
GO:0051607  0.529   2   0.0993700000    0.8377959184    1   12  18  defense response to virus   loc111110197
```

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
