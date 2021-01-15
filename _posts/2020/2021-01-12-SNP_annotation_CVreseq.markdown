---
comments: true
title: CVreseq outlier annotation using snpEff and Annovar
date: '2021-01-12 12:00'
tags:
  - CVreseq
  - ouliter
  - annotation
  - SnpEff
categories:
  - WGS data analysis
---

### SnpEff

SnpEff is a genetic variant annotation and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes). Typically, it answers questions like are SNPs in a gene? In an exon? Do they change protein coding? Do they cause premature stop codons?

Features:

1) Supports over 38,000 genomes       
2) Standard ANN annotation format      
3) Cancer variants analysis      
4) GATK compatible (-o gatk)      
5) HGVS notation      
6) Sequence Ontology standardized terms      

### Reference genome configuration for annotation 

Unfortunately, it does not include the *Crassostrea virginica* genome. Therefore, I have to manually configure the *C. virginica* genome for my own usage. The detailed manual is [here](https://pcingola.github.io/SnpEff/SnpEff_manual.html#run). Below are the detailed steps for genome configuiration,

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

### Annotate the SNPs

```sh
java -jar snpEff.jar eff mygenome 24.recode.vcf > 24_annotation.vcf
```

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

The annotation results (ANN) are listed [here](https://docs.google.com/spreadsheets/d/19TCyi7zKxK7OC6TFrPtCR4DkcgOBqLe20pV2jQw5Jkk/edit?usp=sharing). I searched gene names for the 1st annotation result (1st ANN). 

Note that in this ANN table more than one annotation could be reported for each SNP, and these annotations are sorted by

1) Putative impact: Effects having higher putative impact are first.    
2) Effect type: Effects assumed to be more deleterious effects first.     
3) Canonical transcript before non-canonical.       
4) Marker genomic coordinates (e.g. genes starting before first).        


--- 

### Annovar

Annovar is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants. Four types of the annotation are implemented in the Annovar: 

1) Gene-based annotation: identify whether SNPs or CNVs cause protein coding changes and the amino acids that are affected.     
2) Region-based annotation: identify variants in specific genomic regions, for example, conserved regions among 44 species, predicted transcription factor binding sites, segmental duplication regions, GWAS hits, database of genomic variants, DNAse I hypersensitivity sites, ENCODE H3K4Me1/H3K4Me3/H3K27Ac/CTCF sites, ChIP-Seq peaks, RNA-Seq peaks, or many other annotations on genomic intervals.        
3) Filter-based annotation: identify variants that are documented in specific databases, for example, whether a variant is reported in dbSNP, what is the allele frequency in the 1000 Genome Project, find intergenic variants with GERP++ score<2 or CADD>10, or many other annotations on specific mutations.       
4) Other functionalities: Retrieve the nucleotide sequence in any user-specific genomic positions in batch, identify a candidate gene list for Mendelian diseases from exome data, and other utilities.   

- software

```sh
ANNOVAR  
│  annotate_variation.pl #main function perl script，functions include download the database, and provide three annotation options
│  coding_change.pl #used to infer the changes in protein
│  convert2annovar.pl #format conversion
│  retrieve_seq_from_fasta.pl 
│  table_annovar.pl #annotate script
│  variants_reduction.pl #filter script
├─example #example data
└─humandb #human reference database
```

### Reference genome configuration for annotation 

We need the target reference genome and the annotation file (gtf or gff3 format) for configuration. 

To do that, we need a tool called gtfToGenePred

```sh
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

# convert gtf to GenePred file
gtfToGenePred -genePredExt ref_C_virginica-3.0_top_level.gtf cv30_refGene.txt

# now I renamed the fasta file for Cv genome
mv sequences.fa cv30.fasta #sequences.fa is the Cv genome from NCBI (including the mtDNA)

# build Annovar annotation .fa file
perl retrieve_seq_from_fasta.pl --format refGene --seqfile cv30.fasta cv30_refGene.txt --out cv30_refGeneMrna.fa
NOTICE: Reading region file cv30_refGene.txt ... Done with 67854 regions from 10 chromosomes
NOTICE: Finished reading 11 sequences from cv30.fasta
NOTICE: Finished reading 11 sequences from cv30.fasta
NOTICE: Finished reading 11 sequences from cv30.fasta
NOTICE: Finished reading 11 sequences from cv30.fasta
NOTICE: Finished reading 11 sequences from cv30.fasta
NOTICE: Finished reading 11 sequences from cv30.fasta
NOTICE: Finished reading 11 sequences from cv30.fasta
NOTICE: Finished reading 11 sequences from cv30.fasta
NOTICE: Finished reading 11 sequences from cv30.fasta
NOTICE: Finished reading 11 sequences from cv30.fasta
NOTICE: Finished writting FASTA for 67854 genomic regions to cv30_refGeneMrna.fa
WARNING: 323 gene regions do not have complete ORF (for example, rna27094NC_035783.1:42604388, rna771NC_035780.1:8169903, rna20375NC_035782.1:62198604, rna1433NC_035780.1:13287469, rna22265NC_035783.1:1908843)
```

Now we have two files ready for annotation

```sh
ANNOVAR  
│  cv30_refGeneMrna.fa
│  cv30_refGene.txt
```

### Vcf file formatting

```sh
perl convert2annovar.pl -format vcf4 95.outlier.SNPs.inversion.recode.vcf > 95.outlier.SNPs.inversion.avinput
NOTICE: Finished reading 2365 lines from VCF file
NOTICE: A total of 2269 locus in VCF file passed QC threshold, representing 2269 SNPs (1175 transitions and 1094 transversions) and 0 indels/substitutions
NOTICE: Finished writing 40 SNP genotypes (19 transitions and 21 transversions) and 0 indels/substitutions for 1 sample
```

Code above can only produce 40 SNPs for downstream annotation. Why?

From the Annovar site [https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/input/](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/input/), the author mentioned that "By default, only the first sample in VCF file will be written to output file. The input VCF file contains seven loci, but many of them do not have non-reference genotypes for the first sample, and that is why the output contains only 3 variants." It seems convert2annovar.pl can only format the loci with non-reference genotypes. 

However, with the argument of *-withfreq*, this script could output all SNPs in the original input Vcf file.

```sh
perl convert2annovar.pl -format vcf4 95.outlier.SNPs.inversion.recode.vcf -outfile 95.outlier.SNPs.inversion.avinput -allsample -withfreq
NOTICE: Finished reading 2365 lines from VCF file
NOTICE: A total of 2269 locus in VCF file passed QC threshold, representing 2269 SNPs (1175 transitions and 1094 transversions) and 0 indels/substitutions
NOTICE: Finished writing allele frequencies based on 2269 SNP genotypes (1175 transitions and 1094 transversions) and 0 indels/substitutions for 1 samples
```

Great! This is what I want. Now proceed with the annotation step.

### Annotate the SNPs

```sh
perl annotate_variation.pl -geneanno -dbtype refGene -out 95.outlier.SNPs.inversion -build cv30 95.outlier.SNPs.inversion.avinput ./
# -geneanno: gene based annotation
# -dbtype refGene:  based on refGene database 
# -out zunla: extension file name

perl annotate_variation.pl -geneanno -dbtype refGene -out 95.outlier.SNPs.inversion -build cv30 95.outlier.SNPs.inversion.avinput ./
NOTICE: Output files are written to 95.outlier.SNPs.inversion.variant_function, 95.outlier.SNPs.inversion.exonic_variant_function
NOTICE: Reading gene annotation from ./cv30_refGene.txt ... Done with 67854 transcripts (including 7653 without coding sequence annotation) for 39493 unique genes
NOTICE: Processing next batch with 40 unique variants in 40 input lines
NOTICE: Reading FASTA sequences from ./cv30_refGeneMrna.fa ... Done with 1 sequences
WARNING: A total of 323 sequences will be ignored due to lack of correct ORF annotation
```

### Interpret the annotation results

Two output files are generated: 95.outlier.SNPs.inversion.variant_function and 95.outlier.SNPs.inversion.exonic_variant_function

.variant_function file: the first and second column annotate *variant effects on gene structure and the genes that are affected*, yet the other columns are reproduced from input file

The first column tells whether the variant hit exons or hit intergenic regions, or hit introns, or hit a non-coding RNA genes. If the variant is exonic/intronic/ncRNA, the second column gives the gene name (if multiple genes are hit, comma will be added between gene names); if not, the second column will give the two neighboring genes and the distance to these neighboring genes. [https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/gene/](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/gene/). 

More detailed explanation of variant annotations are given below,

| Value      | Default precedence | Explanation                                                                                                           | Sequence Ontology                          |
|------------|--------------------|-----------------------------------------------------------------------------------------------------------------------|--------------------------------------------|
| exonic     | 1                  | variant overlaps a coding                                                                                             | exon_variant (SO:0001791)                  |
| splicing   | 1                  | variant is within 2-bp of a splicing junction (use -splicing_threshold to change this)                                | splicing_variant (SO:0001568)              |
| ncRNA      | 2                  | variant overlaps a transcript without coding annotation in the gene definition (see Notes below for more explanation) | non_coding_transcript_variant (SO:0001619) |
| UTR5       | 3                  | variant overlaps a 5' untranslated region                                                                             | 5_prime_UTR_variant (SO:0001623)           |
| UTR3       | 3                  | variant overlaps a 3' untranslated region                                                                             | 3_prime_UTR_variant (SO:0001624)           |
| intronic   | 4                  | variant overlaps an intron                                                                                            | intron_variant (SO:0001627)                |
| upstream   | 5                  | variant overlaps 1-kb region upstream of transcription start site                                                     | upstream_gene_variant (SO:0001631)         |
| downstream | 5                  | variant overlaps 1-kb region downtream of transcription end site (use -neargene to change this)                       | downstream_gene_variant (SO:0001632)       |
| intergenic | 6                  | variant is in intergenic region                                                                                       | intergenic_variant (SO:0001628)            |

```sh
head -n 5  95.outlier.SNPs.inversion.variant_function
intronic  gene10029 NC_035782.1 33693099  33693099  A T 0 16827.2 22
intronic  gene10029 NC_035782.1 33693296  33693296  A G 0 11462.4 21
intronic  gene10029 NC_035782.1 33693459  33693459  T G 0 17380.4 31
intronic  gene10029 NC_035782.1 33693556  33693556  A G 0 18352.1 24
intronic  gene10029 NC_035782.1 33693572  33693572  C T 0 19123.6 20
```

.exonic_variant_function file contains the amino acid changes as a result of the exonic variant. Note that only exonic variants are annotated in this file, so the first column gives the line # in the original input file. The second field tells the functional consequences of the variant (possible values in this fields include: nonsynonymous SNV, synonymous SNV, frameshift insertion, frameshift deletion, nonframeshift insertion, nonframeshift deletion, frameshift block substitution, nonframshift block substitution). The third column contains the gene name, the transcript identifier and the sequence change in the corresponding transcript.

More detailed explanation of these exonic_variant_functoin annotations are given below

| Annotation                       | Precedence | Explanation                                                                                                                                                                                                                                                                                         | Sequence Ontology                  |
|----------------------------------|------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------|
| frameshift insertion             | 1          | an insertion of one or more nucleotides that cause frameshift changes in protein coding sequence                                                                                                                                                                                                    | frameshift_elongation (SO:0001909) |
| frameshift deletion              | 2          | a deletion of one or more nucleotides that cause frameshift changes in protein coding sequence                                                                                                                                                                                                      | frameshift_truncation (SO:0001910) |
| frameshift block substitution    | 3          | a block substitution of one or more nucleotides that cause frameshift changes in protein coding sequence                                                                                                                                                                                            | frameshift_variant (SO:0001589)    |
| stopgain                         | 4          | a nonsynonymous SNV, frameshift insertion/deletion, nonframeshift insertion/deletion or block substitution that lead to the immediate creation of stop codon at the variant site. For frameshift mutations, the creation of stop codon downstream of the variant will not be counted as "stopgain"! | stop_gained (SO:0001587)           |
| stoploss                         | 5          | a nonsynonymous SNV, frameshift insertion/deletion, nonframeshift insertion/deletion or block substitution that lead to the immediate elimination of stop codon at the variant site                                                                                                                 | stop_lost (SO:0001578)             |
| nonframeshift insertion          | 6          | an insertion of 3 or multiples of 3 nucleotides that do not cause frameshift changes in protein coding sequence                                                                                                                                                                                     | inframe_insertion (SO:0001821)     |
| nonframeshift deletion           | 7          | a deletion of 3 or mutliples of 3 nucleotides that do not cause frameshift changes in protein coding sequence                                                                                                                                                                                       | inframe_deletion (SO:0001822)      |
| nonframeshift block substitution | 8          | a block substitution of one or more nucleotides that do not cause frameshift changes in protein coding sequence                                                                                                                                                                                     | inframe_variant (SO:0001650)       |
| nonsynonymous SNV                | 9          | a single nucleotide change that cause an amino acid change                                                                                                                                                                                                                                          | missense_variant (SO:0001583)      |
| synonymous SNV                   | 10         | a single nucleotide change that does not cause an amino acid change                                                                                                                                                                                                                                 | synonymous_variant (SO:0001819)    |
| unknown                          | 11         | unknown function (due to various errors in the gene structure definition in the database file)                                                                                                                                                                                                      | sequence_variant (SO:0001060)      |

```sh
head -n 5 95.outlier.SNPs.inversion.exonic_variant_function
line60  nonsynonymous SNV gene10032:rna16980:exon3:c.A434C:p.E145A, NC_035782.1 33735790  33735790  A C 0 20275.3 11
line397 nonsynonymous SNV gene10035:rna16990:exon2:c.A88T:p.I30F, NC_035782.1 34055840  34055840  A T 0 13493.1 23
line415 nonsynonymous SNV gene10035:rna16990:exon3:c.C242T:p.A81V,  NC_035782.1 34075009  34075009  C T 0 15876 20
line511 nonsynonymous SNV gene10035:rna16990:exon6:c.G677A:p.G226D, NC_035782.1 34109058  34109058  G A 0 20133.2 14
line671 nonsynonymous SNV gene10035:rna16990:exon14:c.A1952T:p.N651I, NC_035782.1 34174985  34174985  A T 0 17720 23
```

- Useful options:

1) -csvout  output the csv file (comma-delimited file).      
2) -operation followed by some arguments: g = gene-based; gx = gene-based with cross-reference (need provide -xref argument); r = region-based; f = filter-based.      
3) -withfreq  This is a very important argument. with this option the output file should contain all loci from the input file   
4) 
