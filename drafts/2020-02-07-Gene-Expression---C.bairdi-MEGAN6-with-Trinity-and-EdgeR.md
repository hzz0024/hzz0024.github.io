---
layout: post
title: Gene Expression - C.bairdi MEGAN6 with Trinity and EdgeR
date: '2020-02-07 05:45'
tags:
  - Trinity
  - EdgeR
  - Tanner crab
  - Chionoecetes bairdi
  - gene expression
categories:
  - Miscellaneous
---



---

#### RESULTS

Output folder:

- [20200207_cbai_DEG/](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/)

Comparisons:

---

D12_infected-vs-D12_uninfected

Took a little less than 20mins to run:

![Mox runtime for D12 infected vs D12 uninfeced](https://github.com/RobertsLab/sams-notebook/blob/master/images/screencaps/20200207_cbai_DEG_D12_infected-vs-D12_uninfected_runtime.png?raw=true)

- [D12_infected-vs-D12_uninfected/](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/D12_infected-vs-D12_uninfected)

Only a single DEG, which is upregulated in the infected set:

![MA/volcano plot of D12 infected vs D12 uninfeced](https://github.com/RobertsLab/sams-notebook/blob/master/images/screencaps/20200207_cbai_DEG_D12_infected-vs-D12_uninfected_MA-plot.png?raw=true)

- [salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.DE.subset](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/D12_infected-vs-D12_uninfected/edgeR.24484.dir/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.DE.subset)

TRINITY_DN10191_c0_g1 - [SPID: Q36421](https://www.uniprot.org/uniprot/Q36421) (Cyctochrome c oxidase I)


---

D12_infected-vs-D26_infected

Took ~18mins to run:

![D12 infected vs D26 infected runtime](https://github.com/RobertsLab/sams-notebook/blob/master/images/screencaps/20200207_cbai_DEG_D12_infected-vs-D26_infected_runtime.png?raw=true)

- [D12_infected-vs-D26_infected/](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/D12_infected-vs-D26_infected)

No differentially expressed genes between these two groups.

NOTE: Since no DEGs, that's why this run shows as `FAILED` in the above runtime screencap. This log file captures the error message that kills the job and generates the `FAILED` indicator:

[20200207_cbai_DEG/D12_infected-vs-D26_infected/edgeR.21680.dir/diff_expr_stderr.txt](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/D12_infected-vs-D26_infected/edgeR.21680.dir/diff_expr_stderr.txt)

`Error, no differentially expressed transcripts identified at cuttoffs: P:0.05, C:1 at /gscratch/srlab/programs/trinityrnaseq-v2.9.0/Analysis/DifferentialExpression/analyze_diff_expr.pl line 203.`

---

D12_uninfected-vs-D26_uninfected


Took ~18mins to run:

![D12 uninfected vs D26 uninfected runtime](https://github.com/RobertsLab/sams-notebook/blob/master/images/screencaps/20200207_cbai_DEG_D12_uninfected-vs-D26_uninfected_runtime.png?raw=true)


- [D12_uninfected-vs-D26_uninfected/](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/D12_uninfected-vs-D26_uninfected)

No differentially expressed genes between these two groups.

NOTE: Since no DEGs, that's why this run shows as `FAILED` in the above runtime screencap. This log file captures the error message that kills the job and generates the `FAILED` indicator:

[20200207_cbai_DEG/D12_uninfected-vs-D26_uninfected/edgeR.27147.dir/diff_expr_stderr.txt](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/D12_uninfected-vs-D26_uninfected/edgeR.27147.dir/diff_expr_stderr.txt)

`Error, no differentially expressed transcripts identified at cuttoffs: P:0.05, C:1 at /gscratch/srlab/programs/trinityrnaseq-v2.9.0/Analysis/DifferentialExpression/analyze_diff_expr.pl line 203.`

---


D12-vs-D26

Took ~40mins to run:

![D12 vs D26 runtime](https://github.com/RobertsLab/sams-notebook/blob/master/images/screencaps/20200207_cbai_DEG_D12-vs-D26_runtime.png?raw=true)

- [D12-vs-D26/](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/D12-vs-D26)


---

D26_infected-vs-D26_uninfected

- [D26_infected-vs-D26_uninfected/](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/D26_infected-vs-D26_uninfected)

No differentially expressed genes between these two groups.

NOTE: Since no DEGs, that's why this run shows as `FAILED` in the above runtime screencap. This log file captures the error message that kills the job and generates the `FAILED` indicator:

[20200207_cbai_DEG/D26_infected-vs-D26_uninfected/edgeR.20733.dir/diff_expr_stderr.txt](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/D26_infected-vs-D26_uninfected/edgeR.20733.dir/diff_expr_stderr.txt)

`Error, no differentially expressed transcripts identified at cuttoffs: P:0.05, C:1 at /gscratch/srlab/programs/trinityrnaseq-v2.9.0/Analysis/DifferentialExpression/analyze_diff_expr.pl line 203.`

---

infected-vs-uninfected

Took ~40mins to run:

![infected vs uninfected runtim](https://github.com/RobertsLab/sams-notebook/blob/master/images/screencaps/20200207_cbai_DEG_infected-vs-uninfected_runtime.png?raw=true)

- [infected-vs-uninfected/](https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/infected-vs-uninfected)
