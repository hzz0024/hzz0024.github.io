---
comments: true
title: CVreseq outlier annotation using snpEff and Annovar
date: '2021-01-21 12:00'
tags:
  - CVreseq
  - ouliter
  - annotation
  - SnpEff
categories:
  - WGS data analysis
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

### How to interpret the annotation results

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

### Results

- General annotation for all SNPs in the dataset. This tells whether the variant hit exons or hit intergenic regions, or hit introns, or hit a non-coding RNA genes.

<img src="https://hzz0024.github.io/images/CVreseq_annovar/95.outlier.SNPs.inversion.variant_function.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/CVreseq_annovar/95.outlier.SNPs.no_inversion.variant_function.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/CVreseq_annovar/95.outlier.SNPs.wildae.no_inversion.variant_function.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/CVreseq_annovar/CS_HC-HCVA_CLP.outliers.variant_function.jpg" alt="img" width="800"/>


- Exonic SNPs in each dataset. This tells the functional consequences of the variant.

<img src="https://hzz0024.github.io/images/CVreseq_annovar/95.outlier.SNPs.inversion.exonic_variant_function.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/CVreseq_annovar/95.outlier.SNPs.no_inversion.exonic_variant_function.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/CVreseq_annovar/95.outlier.SNPs.wildae.no_inversion.exonic_variant_function.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/CVreseq_annovar/CS_HC-HCVA_CLP.outliers.exonic_variant_function.jpg" alt="img" width="800"/>

### A comparison between outlier SNPs and its corresponding reference datasets

`Total inversion segments:` 1573 (based on CNV and big inversions)

- Brief summary for number of SNPs

|                                                           | Number of SNPs |
|-----------------------------------------------------------|----------------|
| Genome-wide SNPs within inversions                        | 1125493        |
| Genome-wide SNPs outside inversions                       | 4640654        |
| SNPs in chromosome 2 and outside the inversions           | 668282         |

- Summery of exonic SNP annotation

| Number of SNPs                                                  | Synonymous | Nonsynonymous | Stopgain | Stoploss | Unknown | Total  |
|-----------------------------------------------------------------|------------|---------------|----------|----------|---------|--------|
| Exonic outlier SNPs within inversions                           | 18         | 14            | 0        | 0        | 0       | 32     |
| Exonic outlier SNPs outside inversions                          | 203        | 159           | 0        | 0        | 3       | 365    |
| Exonic outlier SNPs in wild Atlantic   contrasts (no inversion) | 23         | 24            | 0        | 0        | 0       | 47     |
| Exonic outlier SNPs in salinity contrast   (no inversion)       | 186        | 30            | 0        | 0        | 0       | 216    |
| Exonic genome-wide SNPs within inversions                       | 84411      | 39763         | 282      | 34       | 3058    | 127548 |
| Exonic genome-wide SNPs outside inversions                      | 380535     | 172795        | 1276     | 157      | 18830   | 573593 |
| Exonic SNPs in chromosome 2 and outside inversions              | 54505      | 23054         | 155      | 4        | 5312    | 83030  |


- Exonic genome-wide SNPs within inversions vs outlier SNPs within inversions

<img src="https://hzz0024.github.io/images/CVreseq_annovar/95.outlier.SNPs.inversion.jpg" alt="img" width="800"/>

- Exonic genome-wide SNPs outside inversions vs outlier SNPs outside inversions 

<img src="https://hzz0024.github.io/images/CVreseq_annovar/95.outlier.SNPs.noinversion.jpg" alt="img" width="800"/> 

- Exonic outlier SNPs in wild Atlantic contrasts (no inversion) vs genome-wide SNPs outside inversions

<img src="https://hzz0024.github.io/images/CVreseq_annovar/95.outlier.SNPs.wildae.no_inversion.jpg" alt="img" width="800"/> 

- Exonic outlier SNPs in salinity contrast (no inversion) vs genome-wide SNPs outside inversions

<img src="https://hzz0024.github.io/images/CVreseq_annovar/CS_HC-HCVA_CLP_ratio.jpg" alt="img" width="800"/>

- Test for significance

`Exonic genome-wide SNPs outside inversions vs outlier SNPs outside inversions`

|                                           | Synonymous | Nonsynonymous | Marginal Row Totals |
|-------------------------------------------|------------|---------------|---------------------|
| Exonic outlier SNPs outside   inversions  | 203        | 159           | 362                 |
| Expected number based genome-wide data    | 249        | 113           | 362                 |
| Marginal Column Totals                    | 452        | 272           | 724                 |

The chi-square statistic is 12.4608. The p-value is .000416. Significant at p < .05.

`Exonic outlier SNPs in salinity contrast (no inversion) vs genome-wide SNPs outside inversions`

|                                                         | Synonymous | Nonsynonymous | Marginal Row Totals |
|---------------------------------------------------------|------------|---------------|---------------------|
| Exonic outlier SNPs in salinity contrast (no inversion) | 186        | 30            | 216                 |
| Expected number based genome-wide data                  | 149        | 67            | 216                 |
| Marginal Column Totals                                  | 335        | 97            | 432                 |

The chi-square statistic is 18.2. The p-value is .00002. Significant at p < .05

