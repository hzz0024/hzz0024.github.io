---
comments: true
title: Perl script to convert vcf to genepop 
date: '2020-03-12 12:00'
tags:
  - VCF
  - genepop
  - perl
  - script
categories:
  - tools
---

I customized an existing perl script to help convert the vcf to genepop format

The original perl script is here https://github.com/z0on/2bRAD_denovo/blob/master/vcf2genepop.pl

The customized one can be downloaded <a href="https://hzz0024.github.io/scripts/vcf2genepop_hg.pl" download>or here</a>

The only chnage part is to convert allele A, C, G, T to 01, 02, 03, 04

### usage

Arguments: 
vcf=[filename] : vcf file to convert
   pops=[list] : comma-separated perl patterns to determine population 
                 affiliation of samples in the VCF file based on sample names

### Data preparation

vcf files

### Example run

```sh
vcf2genepop.pl vcf=filtered.vcf pops=O,K,M,S > filtered.gen
```