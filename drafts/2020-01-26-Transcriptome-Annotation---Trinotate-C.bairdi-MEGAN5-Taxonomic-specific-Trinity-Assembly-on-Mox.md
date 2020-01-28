---
layout: post
title: Transcriptome Annotation - Trinotate C.bairdi MEGAN5 Taxonomic-specific Trinity Assembly on Mox
date: '2020-01-26 11:16'
tags:
  - Trinotate
  - annotation
  - transcriptome
  - Tanner crab
  - Chionoecetes bairdi
categories:
  - Miscellaneous
---
After performing [_de novo_ assembly on our Tanner crab MEGAN6 taxonomic-specific RNAseq data on 20200122](https://robertslab.github.io/sams-notebook/2020/01/22/Transcriptome-Assembly-C.bairdi-with-MEGAN6-Taxonomy-specific-Reads-with-Trinity-on-Mox.html) and performing [BLASTx annotation on 20200123](https://robertslab.github.io/sams-notebook/2020/01/23/Transcriptome-Annotation-C.bairdi-MEGAN-Trinity-Assembly-Using-DIAMOND-BLASTx-on-Mox.html), I continued the annotation process by running [Trinotate](https://github.com/Trinotate/Trinotate.github.io/wiki).

Trinotate will perform functional annotation of the transcriptome assembly, including GO terms and an annotation feature map that can be used in subsequent Trinity-based differential gene expression analysis so that functional annotations are carried downstream through that process.

SBATCH script (GitHub):

- [20200126_cbai_trinotate_megan.sh](https://github.com/RobertsLab/sams-notebook/blob/master/sbatch_scripts/20200126_cbai_trinotate_megan.sh)



---

#### RESULTS

Run time was ~30mins:

![Cbai Trinotate runtime](https://github.com/RobertsLab/sams-notebook/blob/master/images/screencaps/20200126_cbai_trinotate_megan_runtime.png?raw=true)

Output folder:

- []()
