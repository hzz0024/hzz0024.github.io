---
layout: post
title: RNAseq Reads Extractions - C.bairdi Taxonomic Reads Extractions with MEGAN6 on swoose
date: '2020-01-14 10:36'
tags:
  - Tanner crab
  - MEGAN6
  - taxonomy
  - Chionoecetes bairdi
  - swoose
categories:
  - Miscellaneous
---
[I previously ran BLASTx and "meganized" the output DAA files on 20200103]() and now need to use MEGAN6 to bin the results into the proper taxonomies. This is accomplished using the MEGAN6 graphical user interface (GUI). This is how the process goes:

1. File > Import from BLAST...

2. Select all "meganized" DAA files for a given set of sequencing (e.g. `304428_S1_L001_R1_001.blastx.daa`, `304428_S1_L001_R2_001.blastx.daa`, `304428_S1_L002_R1_001.blastx.daa`, `304428_S1_L001_R2_001.blastx.daa` )

3. Check the "Paired Reads" box. (I don't think this actually does anything, though...)

4. Click "Next"

5. Check the "Analyze Taxonomy Content" box.

6. Click "Load Accession mapping file" and find the mapping file used for "meganizing the DAA file": `prot_acc2tax-Jul2019X1.abin`

7. Click "Next"

8. Click "Load Accession mapping file" and find the mapping file used for "meganizing the DAA file": `acc2eggnog-Jul2019X.abin`

9. Click "Next"

10. Click "Load Accession mapping file" and find the mapping file used for "meganizing the DAA file": `acc2interpro-Jul2019X.abin`

11. Click "Next" twice, to advance to the "SEED" tab.

12. Click "Load Accession mapping file" and find the mapping file used for "meganizing the DAA file": `acc2seed-May2015XX.abin`

13. Click "Next" twice, to advance to the "LCA Params" tab.

14. Click "Apply"

This will initiate the import process and will create a special MEGAN file: RMA6.

NOTE: This will take a long time _and_ will require a significant amount of disk space! The final files aren't ridiculously large, but the intermediate file that gets generated quickly becomes extremely large (i.e. hundreds of GB)!




---

#### RESULTS

Output folder:

- []()
