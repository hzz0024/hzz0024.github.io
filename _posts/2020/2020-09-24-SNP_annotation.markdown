---
comments: true
title: DelBay19 outlier annotation using snpEff
date: '2020-09-24 12:00'
tags:
  - DelBay19
  - ouliter
  - annotation
  - WGS
  - Fisher's exact
  - 
categories:
  - WGS data analysis
---

After identifing some potential outliers using Fisher's approach, I conducted marker annotation using a program called snpeff. 

snpeff is a genetic variant annotation and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).

Features:

1) Supports over 38,000 genomes       
2) Standard ANN annotation format      
3) Cancer variants analysis      
4) GATK compatible (-o gatk)      
5) HGVS notation      
6) Sequence Ontology standardized terms      

---

### Reference genome configuration for annotation 

Unfortunately, it does not include the Crassostrea virginica genome. Therefore, I have to manually configure the C. virginica genome for my own usage. The detailed manual is [here](https://pcingola.github.io/SnpEff/SnpEff_manual.html#run). Below are the detailed steps for genome configuiration,

1) Add a genome to the configuration file

```sh
vi snpEffect.config 
# Add the following lines below the Databases & Genomes section
#-------------------------------------------------------------------------------
# Databases & Genomes
#
# One entry per genome version.
#
# For genome version 'ZZZ' the entries look like
#   ZZZ.genome              : Real name for ZZZ (e.g. 'Human')
#   ZZZ.reference           : [Optional] Comma separated list of URL to site/s where information for building ZZZ database was extracted.
#   ZZZ.chrName.codonTable  : [Optional] Define codon table used for chromosome 'chrName' (Default: 'codon.Standard')
#
#-------------------------------------------------------------------------------

# my genome
mygenome.genome : Crassostrea_virginica
```
2) Building a database from the GTF files

```sh
# rename the gtf file 
mv ref_C_virginica-3.0_nDNA.gtf genes.gtf

# rename the fasta file
mv cv30.fa sequences.fa

# create the genome folder in the data 
mkdir mygenome
cp genes.gtf data/mygenome
cp sequences.fa data/mygenome

Note: the gtf and fa files must be named as genes.gtf and sequences.fa
```

3) Using snpEff.jar to build the refernece genome

```sh
java -jar snpEff.jar build -gtf22 -v mygenome
```

---

### Vcf file prepariation

- Shell script to convert the chromosome id

```sh
for i in *.vcf; do
sed -i.bak 's/NC_035780.1/1/g;s/NC_035781.1/2/g;s/NC_035782.1/3/g;s/NC_035783.1/4/g;s/NC_035784.1/5/g;s/NC_035785.1/6/g;s/NC_035786.1/7/g;s/NC_035787.1/8/g;s/NC_035788.1/9/g;s/NC_035789.1/10/g;s/NC_007175.2/11/g' $i
done
```

- Python script useful to add the SNP id to the vcf file (produced from Angsd)

```python
# change the fname to target vcf file
fname = 'XXXXXXXX.vcf'
outname = fname + '.out'

idx = 0
with open(fname, 'r') as f, open(outname, 'w') as w:
    for l in f:
        if l.startswith('#'):
            pass
        else:
            idx += 1
            ss = l.split()
            chrom = ss[0]
            pos = ss[1]
            ID = chrom + '_' + pos
            assert ss[2] == '.'
            ss[2] = ID
            l = '\t'.join(ss)
            l += '\n'
        w.write(l)
```

- Using vcftools to extract the target SNPs

```sh
vcftools --vcf ALL_maf0.05_pctind0.7_cv30_reformat.vcf --snps Fisher_fdr_0.05_24_snp.list --recode --recode-INFO-all --out 24
```

---

### Annotate the SNPs

```sh
java -jar snpEff.jar eff mygenome 24.recode.vcf > 24_annotation.vcf
```

--- 

### Interpret the annotation results

|Chromosome	 | Length    | Variants |  Variants rate|
|------------|-----------|----------|---------------|
|1	         |65,668,440 |	6	    | 10,944,740    |
|3	         |77,061,148 |	2	    | 38,530,574    |
|4	         |59,691,872 |	2	    | 29,845,936    |
|5        	 |98,698,416 |	10      | 9,869,841     |
|7	         |57,830,854 |	1	    | 57,830,854    |
|8	         |75,944,018 |	2	    | 37,972,009    |
|10	         |32,650,045 |	1	    | 32,650,045    |
|Total	     |467,544,793|	24	    | 19,481,033    |

Number of effects: 96. Here Effects are categorized by 'impact': {High, Moderate, Low, Modifier}. A list of "impact" are pre-defined categories to help users find more significant variants (see [snpeff input/output files](https://pcingola.github.io/SnpEff/se_inputoutput/#effect-prediction-details)). Below are the tables for "Number of effects" and "Number of effects by region"

<img src="https://hzz0024.github.io/images/SNP_annotation/impact.jpg" alt="img" width="800"/>

- Number of effects

|Type (alphabetical order)	|	Count |	Percent |
|---------------------------|---------|---------|
|LOW	                    |  	7	  | 7.292%  |
|MODERATE	                |	4	  | 4.167%  |
|MODIFIER	 	            |  85	  |88.542%  |

- Number of effects by region (table followed by figure)

|Type (alphabetical order)	| 	Count |	Percent |
|---------------------------|---------|---------|
|DOWNSTREAM	                |	13	  | 13.542% |
|EXON	 	                |   12	  |  12.5%  |
|INTERGENIC	 	            |   5	  |  5.208% |
|INTRON	 	                |   49	  | 51.042% |
|UPSTREAM	 	            |   12	  | 12.5%   |
|UTR_3_PRIME	            | 	4	  | 4.167%  |
|UTR_5_PRIME	 	        |   1	  | 1.042%  |

<img src="https://hzz0024.github.io/images/SNP_annotation/DelBay19_Fish_fdr0.05_outliers.jpg" alt="img" width="800"/>

Here I only focused on the mumber of effects in the exon regions, which is 12. 

|Type (alphabetical order) |   Count  | Percent |
|--------------------------|----------|---------|
|MISSENSE                  |    4     |	33.33%  |
|NONSENSE                  |    1     |	 8.33%  |
|SILENT	                   |    7     | 58.33%  |

The annotation results (ANN) are listed [here](https://docs.google.com/spreadsheets/d/19TCyi7zKxK7OC6TFrPtCR4DkcgOBqLe20pV2jQw5Jkk/edit?usp=sharing)

Note that in the table more than one annotation could be reported for each SNP, and these annotations are sorted by

1) Putative impact: Effects having higher putative impact are first.    
2) Effect type: Effects assumed to be more deleterious effects first.     
3) Canonical transcript before non-canonical.       
4) Marker genomic coordinates (e.g. genes starting before first).        

