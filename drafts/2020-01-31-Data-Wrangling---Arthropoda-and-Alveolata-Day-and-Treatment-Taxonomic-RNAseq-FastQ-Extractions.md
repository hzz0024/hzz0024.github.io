---
layout: post
title: Data Wrangling - Arthropoda and Alveolata Day and Treatment Taxonomic RNAseq FastQ Extractions
date: '2020-01-31 08:33'
tags:
  - hematodinium
  - Tanner crab
  - Chionoecetes bairdi
  - MEGAN6
categories:
  - Miscellaneous
---
After using MEGAN6 to [extract _Arthropoda_ and _Alveolata_ reads from our RNAseq data on 20200114](https://robertslab.github.io/sams-notebook/2020/01/14/RNAseq-Reads-Extractions-C.bairdi-Taxonomic-Reads-Extractions-with-MEGAN6-on-swoose.html), I had then extracted taxonomic-specific reads and aggregated each into basic Read 1 and Read 2 FastQs to simplify [transcriptome assembly for _C.bairdi_](https://robertslab.github.io/sams-notebook/2020/01/22/Transcriptome-Assembly-C.bairdi-with-MEGAN6-Taxonomy-specific-Reads-with-Trinity-on-Mox.html) and [for _Hematodinium_](https://robertslab.github.io/sams-notebook/2020/01/22/Transcriptome-Assembly-Hematodinium-with-MEGAN6-Taxonomy-specific-Reads-with-Trinity-on-Mox.html). That was fine and all, but wasn't fully thought through.

For gene expression analysis, I need the FastQs based on infection status and sample days. So, I need to modify the read extraction procedure to parse reads based on those conditions. I could've/should've done this originally, as I could've just assembled the transcriptome from the FastQs I'm going to generate now. Oh well.

Anyway, here's a brief rundown of the approach:

1. Create list of unique read headers from MEGAN6 FastA files.

2. Use list with `seqtk` program to pull out corresponding FastQ reads from the trimmed FastQ R1 and R2 files.

The entire procedure is documented in a Jupyter Notebook below.

Jupyter notebook (GitHub):

- [20200131_swoose_cbai_megan_day-treatment_read_extractions.ipynb](https://github.com/RobertsLab/code/blob/master/notebooks/sam/20200131_swoose_cbai_megan_day-treatment_read_extractions.ipynb)

---

#### RESULTS

Output folders:

- [20200131.C_bairdi_megan_reads](https://gannet.fish.washington.edu/Atumefaciens/20200131.C_bairdi_megan_reads/)

- [20200131.Hematodinium_megan_reads/](https://gannet.fish.washington.edu/Atumefaciens/20200131.Hematodinium_megan_reads/)


We now have to distinct sets of RNAseq reads from _C.bairdi_ (_Arhtropoda_) and _Hematodinium_ (_Alveolata_), split by infection status and sample day! Will get some gene expression analysis going.
