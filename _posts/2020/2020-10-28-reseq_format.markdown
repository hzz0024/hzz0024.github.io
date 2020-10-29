---
comments: true
title: Resequencing data formating and examination
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

- Check common snps between thinned and masked vcf files

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
```

3) extract the common snps in chromosome 2 for domestic-wild populations 

```sh
grep "2_" common.txt | wc -l
29220
grep "2_" common.txt > common_chr2.txt # snp list for next step

vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf --keep SL --snps common_chr2.txt --chr 2 --recode --recode-INFO-all --out SL_chr2_mask

vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf --keep OBOYS2 --snps common_chr2.txt --chr 2 --recode --recode-INFO-all --out OBOYS2_chr2_mask

vcftools --vcf SNP.ORIGINAL.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf --keep SL --snps common_chr2.txt --chr 2 --recode --recode-INFO-all --out SL_chr2_original

vcftools --vcf SNP.ORIGINAL.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf --keep OBOYS2 --snps common_chr2.txt --chr 2 --recode --recode-INFO-all --out OBOYS2_chr2_original
```

- Estimate Pi and Tajima's D

```sh
#Pi 
# for single snp diversity
for pop in OBOYS2_chr2_mask SL_chr2_mask OBOYS2_chr2_original SL_chr2_original; do
        vcftools --vcf $pop'.recode.vcf' --site-pi --out $pop
done
# diversity in windows
for pop in OBOYS2_chr2_mask SL_chr2_mask OBOYS2_chr2_original SL_chr2_original; do
    for win in 100 200 500 1000 5000; do
        vcftools --vcf $pop'.recode.vcf' --window-pi $win --out $pop'_'$win
    done
done

#Tajima'D
# bins in specific size
for pop in OBOYS2_chr2_mask SL_chr2_mask OBOYS2_chr2_original SL_chr2_original; do
    for win in 100 200 500 1000 5000; do
        vcftools --vcf $pop'.recode.vcf' --TajimaD $win --out $pop'_D'$win
    done
done
```

---

Now let us move back to our questions:

#### Question 1: how many SNPs in Katie’s LD-thinned SNP list are also SNPs in the “masked” vcfs? If most of them are still SNPs in “masked” then we can do some preliminary comparisons with the intersection SNPs that occur in both.

LD-thinned SNP list: 334011      
Common SNPs in the vcf with masked genome: 240908      
Percent: 240908/334011 = 72.13%     

#### Question 2: It would be great to see, for a single chromosome with strong haplotig content (LG 1-5), a plot of sliding window pie and Tajima’s D, compared between original and masked datasets (compared for one wild population and also for one related domesticated strain).

Target chromosome: chromosome 2   
Target domestic vs wild contrast: OBOYS2 (Louisiana selected line) vs. SL (Louisiana wild line)

- pi compared between original and masked datasets

<img src="https://hzz0024.github.io/images/CVreseq/SL_chr2_5000.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/CVreseq/OBOYS2_chr2_5000.jpg" alt="img" width="800"/>




