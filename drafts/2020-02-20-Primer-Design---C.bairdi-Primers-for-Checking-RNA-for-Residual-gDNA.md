---
layout: post
title: Primer Design - C.bairdi Primers for Checking RNA for Residual gDNA
date: '2020-02-20 15:01'
tags:
  - Primer3
  - Tanner crab
  - Chionoecetes bairdi
  - transcriptome
  - DEGs
  - Trinity
  - Trinotate
categories:
  - Miscellaneous
---
Getting ready to run some qPCRs and first we need to confirm that our RNA is actually DNA-free. Before we can do that, we need some primers to use, so I decided to semi-arbitrarily select three different gene targets from our [MEGAN6 taxonomic-specific Trinity assembly from 20200122](https://robertslab.github.io/sams-notebook/2020/01/22/Transcriptome-Assembly-C.bairdi-with-MEGAN6-Taxonomy-specific-Reads-with-Trinity-on-Mox.html).

I used our [recent differential gene expression analysis]() to identify those genes which were highly differentially expressed in infected vs. uninfected samples.

Overall, the process went something like this:



---

#### RESULTS

```
PRIMER PICKING RESULTS FOR cbai_TRINITY_DN6411_c0_g2_i1

No mispriming library specified
Using 1-based sequence positions
OLIGO            start  len      tm     gc%   any    3' seq
LEFT PRIMER         80   20   59.94   45.00  4.00  0.00 TGCCGGTAAGGTGAAAAATC
RIGHT PRIMER       261   20   59.97   45.00  2.00  2.00 AAATCCGCAACCAATACAGC
SEQUENCE SIZE: 334
INCLUDED REGION SIZE: 334

PRODUCT SIZE: 182, PAIR ANY COMPL: 3.00, PAIR 3' COMPL: 0.00

    1 GTTTTTTCCTTTTTCGTTTTCTACATATATTAACCCCCCTTTATTAAACAATGGGTAAAG


   61 TCCACGGTTCCTTGGCTCGTGCCGGTAAGGTGAAAAATCAGACCCCGAAAGTTGCCAAGA
                         >>>>>>>>>>>>>>>>>>>>                     

  121 TGGAGAAGAAGAAGTCTCTCACGGGCCGCGCCAAGAAACGCATGCAGTACAACCGTCGTT


  181 TCGTGAACATCGTGCGGGCAGGTGGCCCCAAGCGCGGCCCTAATTCCAACCAGAAGTAAA


  241 GGCTGTATTGGTTGCGGATTTTAGGTGTTAACGATGCGCTGGACTTCCTCCTCTATATGA
       <<<<<<<<<<<<<<<<<<<<                                       

  301 GTATCATGGGATGGATGCAACGAACTTGATGGAC
```

---

```
	PRIMER PICKING RESULTS FOR cbai_TRINITY_DN13073_c0_g1_i1

No mispriming library specified
Using 1-based sequence positions
OLIGO            start  len      tm     gc%   any    3' seq
LEFT PRIMER         65   20   60.29   50.00  4.00  0.00 CGAGTGTTTCCAAGCCTGTT
RIGHT PRIMER       215   20   60.07   50.00  4.00  0.00 GTGAATACGCCTTCCTTCCA
SEQUENCE SIZE: 237
INCLUDED REGION SIZE: 237

PRODUCT SIZE: 151, PAIR ANY COMPL: 3.00, PAIR 3' COMPL: 1.00

    1 GTAGTATTCTGGAATCGGCGTTTTTTGTTTGTGTAATCCGTGGAAATGGACATATCTCAA


   61 CCCGCGAGTGTTTCCAAGCCTGTTTTTACACGCTTGACCGACCTCGCGAGCGAACTGCTC
          >>>>>>>>>>>>>>>>>>>>                                    

  121 GGCTCGAAGGTGCTTTTTGCCACCGATCAGTGGTTTGCCGAAGCTTCAAATTTACTCAAG


  181 AGTGAAGAGCCGGTATGGAAGGAAGGCGTATTCACCGAACATGGAAAATGGATGGAC
                     <<<<<<<<<<<<<<<<<<<<                      
```

---

```
PRIMER PICKING RESULTS FOR cbai_TRINITY_DN6549_c0_g1_i1

No mispriming library specified
Using 1-based sequence positions
OLIGO            start  len      tm     gc%   any    3' seq
LEFT PRIMER        297   20   60.00   45.00  4.00  2.00 CGGTTTGTTTGAACGGCTAT
RIGHT PRIMER       577   20   59.95   50.00  4.00  3.00 GATAAAGCTCGGCATTCTGC
SEQUENCE SIZE: 647
INCLUDED REGION SIZE: 647

PRODUCT SIZE: 281, PAIR ANY COMPL: 3.00, PAIR 3' COMPL: 0.00

    1 TGCGGGAATATCTTTAAATACTATATACTCGGGTAGCGTCTTGGAATGTCATGTGAGGGA


   61 AATTCAGACCCGCACCATGATTATCAGGCATCCCTGAACCAGCAAGATGCGATCCGGCAG


  121 GAAGCGTCCGTCGATCACCCGTTGATGAAGAAGCGCGAGCCCGTAGGGGCATCGCTGAAC


  181 GAGCAGTTCGCGGAGAATAAGAACTTCCTACAGAAGGTCGCTTCAATCGCGGCCAAGTAT


  241 GAGTTCATTCGACGGGCGAGACCGGACGGCAATTGCTTTTACCGCACGTATCTGTTCGGT
                                                              >>>>

  301 TTGTTTGAACGGCTATTGGGCATGTCCCGCGAGGAGCGGGACAAATTTGTCGTGTTTCTC
      >>>>>>>>>>>>>>>>                                            

  361 AAGAAATCACTGGATGATGTGCTTTGCCAAGGGTATGAGCGATTTGCGGTAGAAGAAATG


  421 CACGAAGATATCCTTGAAGAGTTTGAGAAACTCGCTCAGAATGACAATGCAACCGTCGGC


  481 GATATCGAGACGATATTCGACGAGGAAAGGCATTACCATATTTGCTACTTGAGGTGCCTA


  541 GCGTCGGCGTACCTCAAGCAGAATGCCGAGCTTTATCAATCGTTCCTCGAAGGCTATGCG
                       <<<<<<<<<<<<<<<<<<<<                       

  601 ACTATAGCAGAGTTCTGCGCTCATGAAGTGGATCCTATGTGGCGCGG
```
