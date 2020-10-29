---
comments: true
title: Resequencing data analysis
date: '2020-10-28 12:00'
tags:
  - Reseq
  - vcf 
categories:
  - WGS data analysis
---

Data overview

- Number of SNPs

```sh
grep -v "^#" SNP.ORIGINAL.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf |wc -l
grep -v "^#" SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf |wc -l 
grep -v "^#" Thinned.SNP.TRSdp5g1FnDNAmaf052alleles.thinnedMatrixAndMetaData5000Window_exclude_LM.vcf |wc -l
```
Original: 4734919 (contracted vcf file size: 13G)
Masked: 6413937 (contracted vcf file size: 17G)
Thinned vcf (old): 334011 

- VCF dataset formating

```sh
1) replace the chromosome ID with numbers 
for i in *.vcf; do
sed -i.bak 's/NC_035780.1/1/g;s/NC_035781.1/2/g;s/NC_035782.1/3/g;s/NC_035783.1/4/g;s/NC_035784.1/5/g;s/NC_035785.1/6/g;s/NC_035786.1/7/g;s/NC_035787.1/8/g;s/NC_035788.1/9/g;s/NC_035789.1/10/g;s/NC_007175.2/11/g' $i
done

2) sort the vcf file based on chromosome and position 
cat SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.sort.vcf

cat SNP.ORIGINAL.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > SNP.ORIGINAL.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.sort.vcf

3) add SNP ID in the vcf file (python script)
fname = 'SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.sort.vcf'
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

fname = 'SNP.ORIGINAL.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.sort.vcf'
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

fname = 'Thinned.SNP.TRSdp5g1FnDNAmaf052alleles.thinnedMatrixAndMetaData5000Window_exclude_LM.sort.vcf'
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

4) rename the vcf files
mv SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.sort.vcf.out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf
mv SNP.ORIGINAL.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.sort.vcf.out SNP.ORIGINAL.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf
mv Thinned.SNP.TRSdp5g1FnDNAmaf052alleles.thinnedMatrixAndMetaData5000Window_exclude_LM.sort.vcf.out Thinned.SNP.TRSdp5g1FnDNAmaf052alleles.thinnedMatrixAndMetaData5000Window_exclude_LM.format.vcf
```

- check common snps between thinned and masked vcf files

1) check common snps using vcftools

```sh
vcftools --gzvcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf --gzdiff Thinned.SNP.TRSdp5g1FnDNAmaf052alleles.thinnedMatrixAndMetaData5000Window_exclude_LM.format.vcf --diff-site --out vcftools_diff_mask_thinned

Found 240908 sites common to both files.
Found 6173029 sites only in main file.
Found 93103 sites only in second file.
Found 0 non-matching overlapping sites.
After filtering, kept 6413937 out of a possible 6413937 Sites
Run Time = 484.00 seconds
```

2) step above will produce a snp list called "vcftools_diff_mask_thinned.diff.sites_in_files", we can extract common snps in this output file (python script)

```py
# extrac the common shared snp
filename = 'vcftools_diff_mask_thinned.diff.sites_in_files'

with open(filename, 'r') as f, open('common.txt', 'w') as w:
    for l in f:
        ss = l.split()
        if ss[3] == 'B': # B stands for Both
            w.write(ss[0]+'_'+ss[1]+'\n')

filename = 'vcftools_diff_mask_thinned.diff.sites_in_files'
# extrac the mask only snps
with open(filename, 'r') as f, open('mask_only.txt', 'w') as w:
    for l in f:
        ss = l.split()
        if ss[3] == '1': # 1 means only show in the masked vcf file
            w.write(ss[0]+'_'+ss[1]+'\n')
# extrac the thinned only snps
with open(filename, 'r') as f, open('thinned_only.txt', 'w') as w:
    for l in f:
        ss = l.split()
        if ss[3] == '2': # 2 means only show in the thinned vcf file
            w.write(ss[0]+'_'+ss[1]+'\n')
```

3) extract the common snps


