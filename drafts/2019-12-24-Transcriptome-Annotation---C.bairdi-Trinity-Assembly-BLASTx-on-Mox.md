---
layout: post
title: Transcriptome Annotation - C.bairdi Trinity Assembly BLASTx on Mox
date: '2019-12-24 20:36'
tags:
  - transcriptome
  - annotation
  - tanner crab
  - blastx
  - mox
  - Chionoecetes bairdi
categories:
  - Miscellaneous
---



---

#### RESULTS

This took a bit over 17hrs to complete:

![cbai blastx runtime screencap](https://github.com/RobertsLab/sams-notebook/blob/master/images/screencaps/20191224_cbai_blastx_outfmt-11_runtime.png?raw=true)

Output folder:

- [20191224_cbai_blastx_outfmt-11/](https://gannet.fish.washington.edu/Atumefaciens/20191224_cbai_blastx_outfmt-11/)

BLASTx output (ASN):

- [20191224_cbai_blastx_outfmt-11/20191224-20191218.C_bairdi.Trinity.fasta.blastx.asn](https://gannet.fish.washington.edu/Atumefaciens/20191224_cbai_blastx_outfmt-11/20191224-20191218.C_bairdi.Trinity.fasta.blastx.asn) (1.9GB)

This BLASTx file will be used for complete transcriptome annotation using Trinotate.

BLAST is also capable of converting the ASN output format into any other BLAST output format (e.g. format 6), so having BLAST results in ASN format provides a bit more flexibility. However, conversion from ASN to another BLAST output format does take a while.
