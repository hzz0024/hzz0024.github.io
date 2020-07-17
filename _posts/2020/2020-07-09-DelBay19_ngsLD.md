---
comments: true
title: DelBay19 ngsLD 
date: '2020-07-17 12:00'
tags:
  - DelBay19
  - deltap
  - ouliter
  - ngsLD
  - WGS
categories:
  - WGS data analysis
---

This post shows the results from ngsLD. The detailed steps for data creation is in post [Rebuild DelBay19 data](https://hzz0024.github.io/2020/07/08/DelBay_data_redo.html)

Briefly, I created ngsLD output for each population: CH, REF, HC, ARN, COH, SR, and NB. For each population, the LD decay curvy will be ploted for chromosome 1-10. 

Initially, I evaluated what is the optimal setting of --max_kb_dist. This means the maximum distance between SNPs (in Kb) to calculate LD. I created three input with 1, 5, and 10 kb as the mximum distance between SNPs for LD calculation.

- REF 1kb

<img src="https://hzz0024.github.io/images/ngsLD/REF_1k.jpg" alt="img" width="800"/>

- REF 5kb

<img src="https://hzz0024.github.io/images/ngsLD/REF_5k.jpg" alt="img" width="800"/>

- REF 10kb

<img src="https://hzz0024.github.io/images/ngsLD/REF_10k.jpg" alt="img" width="800"/>

It looks the curvy became flat after 5kb. I then create the LD decay curvy for each of the examined population with --max_kb_dist 5

- CH

<img src="https://hzz0024.github.io/images/ngsLD/CH_5k.jpg" alt="img" width="800"/>

- REF

<img src="https://hzz0024.github.io/images/ngsLD/REF_5k.jpg" alt="img" width="800"/>

- HC

<img src="https://hzz0024.github.io/images/ngsLD/HC_5k.jpg" alt="img" width="800"/>

- ARN

<img src="https://hzz0024.github.io/images/ngsLD/ARN_5k.jpg" alt="img" width="800"/>

- COH

<img src="https://hzz0024.github.io/images/ngsLD/COH_5k.jpg" alt="img" width="800"/>

- SR

<img src="https://hzz0024.github.io/images/ngsLD/SR_5k.jpg" alt="img" width="800"/>

- NB

<img src="https://hzz0024.github.io/images/ngsLD/NB_5k.jpg" alt="img" width="800"/>

- Estimates with a single population and sample LD at high resolution between 0 and 500 bp

- max 500, sample every 50 bp

```sh
Rscript --vanilla --slave fit_LDdecay.R --ld_files LD.list --out REF_500_25bp.jpg --fit_bin_size 50 --plot_bin_size 50 --max_kb_dist 0.5 --fit_level 3
Random seed: 44970
==> Fitting r2 LD decay assuming a three (rate of decay, max LD and min LD) parameter decay model
Warning message:
Fitting of LD decay is highly unreliable at short distances (<50kb).
      File LD  DecayRate     LDmax     LDmin
1 REF_chr1 r2 0.01824830 0.2629698 0.1129119
2 REF_chr2 r2 0.01767198 0.2706954 0.1088194
```

<img src="https://hzz0024.github.io/images/ngsLD/REF_500_50bp.jpg" alt="img" width="800"/>

- max 300, sample every 25 bp

```sh
Rscript --vanilla --slave fit_LDdecay.R --ld_files LD.list --out REF_300_25bp.jpg --fit_bin_size 25 --plot_bin_size 25 --max_kb_dist 0.3 --fit_level 3
==> Fitting r2 LD decay assuming a three (rate of decay, max LD and min LD) parameter decay model
Warning message:
Fitting of LD decay is highly unreliable at short distances (<50kb).
      File LD  DecayRate     LDmax     LDmin
1 REF_chr1 r2 0.03344789 0.2979377 0.1302768
2 REF_chr2 r2 0.03105276 0.3071178 0.1261241
```

<img src="https://hzz0024.github.io/images/ngsLD/REF_300_25bp.jpg" alt="img" width="800"/>

- max 100, sample every 10 bp

```sh
Rscript --vanilla --slave fit_LDdecay.R --ld_files LD.list --out REF_100_10bp.jpg --fit_bin_size 10 --plot_bin_size 10 --max_kb_dist 0.1 --fit_level 3
==> Fitting r2 LD decay assuming a three (rate of decay, max LD and min LD) parameter decay model
Warning message:
Fitting of LD decay is highly unreliable at short distances (<50kb).
      File LD  DecayRate     LDmax     LDmin
1 REF_chr1 r2 0.06235471 0.3453660 0.1505267
2 REF_chr2 r2 0.05653356 0.3537681 0.1490709
```

<img src="https://hzz0024.github.io/images/ngsLD/REF_100_10bp.jpg" alt="img" width="800"/>

- max 50, sample every 5 bp

```sh
Rscript --vanilla --slave fit_LDdecay.R --ld_files LD.list --out REF_50_5bp.jpg --fit_bin_size 5 --plot_bin_size 5 --max_kb_dist 0.05 --fit_level 3
==> Fitting r2 LD decay assuming a three (rate of decay, max LD and min LD) parameter decay model
Warning message:
Fitting of LD decay is highly unreliable at short distances (<50kb).
      File LD DecayRate     LDmax     LDmin
1 REF_chr1 r2 0.1881567 0.3923010 0.1966648
2 REF_chr2 r2 0.1614777 0.3991607 0.1975584
```

<img src="https://hzz0024.github.io/images/ngsLD/REF_50_5bp.jpg" alt="img" width="800"/>

Note that the LD decay is the averaged restuls, single point r^2 could be as high as ~1, for example, in REF chr 1,

| snp 1           | snp 2              | distance(bp) |   r^2    |
|-----------------|--------------------|--------------|----------|
|NC_035780.1:30466|  NC_035780.1:30473 |      7       | 0.999307 |
|NC_035780.1:30495|  NC_035780.1:30567 |      72      | 0.445740 |
|NC_035780.1:30684|  NC_035780.1:30759 |      75      | 0.515723 |