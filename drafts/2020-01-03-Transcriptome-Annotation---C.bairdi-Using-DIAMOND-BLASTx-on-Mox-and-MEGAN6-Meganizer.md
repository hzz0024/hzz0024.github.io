---
layout: post
title: Transcriptome Annotation - C.bairdi Using DIAMOND BLASTx on Mox and MEGAN6 Meganizer
date: '2020-01-03 10:05'
tags:
  - tanner crab
  - mox
  - MEGAN
  - DIAMOND
  - BLASTx
  - meganizer
  - Chionoecetes bairdi
categories:
  - Miscellaneous
---


SBATCH script (GitHub):

- [20200103_cbai_diamond_blastx.sh](https://github.com/RobertsLab/sams-notebook/blob/master/sbatch_scripts/20200103_cbai_diamond_blastx.sh)

```shell
#!/bin/bash
## Job Name
#SBATCH --job-name=cbai_blastx_DIAMOND
## Allocation Definition
#SBATCH --account=coenv
#SBATCH --partition=coenv
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=20-00:00:00
## Memory per node
#SBATCH --mem=120G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samwhite@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/samwhite/outputs/20200103_cbai_diamond_blastx

## Perform DIAMOND BLASTx on trimmed Chionoecetes bairdi (Tanner crab) FastQ files.

## Trimmed FastQ files originated here:
## https://gannet.fish.washington.edu/Atumefaciens/20191218_cbai_fastp_RNAseq_trimming/

# Exit script if any command fails
set -e

# Load Python Mox module for Python module availability

module load intel-python3_2017

# SegFault fix?
export THREADS_DAEMON_MODEL=1

# Document programs in PATH (primarily for program version ID)

{
date
echo ""
echo "System PATH for $SLURM_JOB_ID"
echo ""
printf "%0.s-" {1..10}
echo "${PATH}" | tr : \\n
} >> system_path.log


# Program paths
diamond=/gscratch/srlab/programs/diamond-0.9.29/diamond

# DIAMOND NCBI nr database
dmnd=/gscratch/srlab/blastdbs/ncbi-nr-20190925/nr.dmnd


# FastQ files directory
fastq_dir=/gscratch/srlab/sam/data/C_bairdi/RNAseq/


# Loop through FastQ files, log filenames to fastq_list.txt.
# Run DIAMOND on each FastQ
for fastq in ${fastq_dir}*fastp-trim*.fq.gz
do
	# Log input FastQs
	echo "${fastq}" >> fastq_list.txt

	# Strip leading path and extensions
	no_path=$(echo "${fastq##*/}")
	no_ext=$(echo "${no_path%%.*}")

	# Run DIAMOND with blastx
	# Output format 100 produces a DAA binary file for use with MEGAN
	${diamond} blastx \
	--db ${dmnd} \
	--query "${fastq}" \
	--out "${no_ext}".blastx.daa \
	--outfmt 100 \
	--top 5 \
	--block-size 15.0 \
	--index-chunks 4
done
```

---

#### RESULTS

Output folder:

- [20200103_cbai_diamond_blastx/](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/)




Here's the full list of output DIAMOND `daa` files and their sizes (note: they're _huge_ files):

- [304428_S1_L001_R1_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/304428_S1_L001_R1_001.blastx.daa) (56GB)

- [304428_S1_L001_R2_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/304428_S1_L001_R2_001.blastx.daa) (54GB)

- [304428_S1_L002_R1_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/304428_S1_L002_R1_001.blastx.daa) (54GB)

- [304428_S1_L002_R2_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/304428_S1_L002_R2_001.blastx.daa) (52GB)

- [329774_S1_L001_R1_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329774_S1_L001_R1_001.blastx.daa) (39GB)

- [329774_S1_L001_R2_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329774_S1_L001_R2_001.blastx.daa) (36GB)

- [329774_S1_L002_R1_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329774_S1_L002_R1_001.blastx.daa) (34GB)

- [329774_S1_L002_R2_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329774_S1_L002_R2_001.blastx.daa) (32GB)

- [329775_S2_L001_R1_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329775_S2_L001_R1_001.blastx.daa) (40GB)

- [329775_S2_L001_R2_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329775_S2_L001_R2_001.blastx.daa) (36GB)

- [329775_S2_L002_R1_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329775_S2_L002_R1_001.blastx.daa) (37GB)

- [329775_S2_L002_R2_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329775_S2_L002_R2_001.blastx.daa) (32GB)

- [329776_S3_L001_R1_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329776_S3_L001_R1_001.blastx.daa) (35GB)

- [329776_S3_L001_R2_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329776_S3_L001_R2_001.blastx.daa) (32GB)

- [329776_S3_L002_R1_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329776_S3_L002_R1_001.blastx.daa) (30GB)

- [329776_S3_L002_R2_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329776_S3_L002_R2_001.blastx.daa) (29GB)

- [329777_S4_L001_R1_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329777_S4_L001_R1_001.blastx.daa) (40GB)

- [329777_S4_L001_R2_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329777_S4_L001_R2_001.blastx.daa) (34GB)

- [329777_S4_L002_R1_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329777_S4_L002_R1_001.blastx.daa) (36GB)

- [329777_S4_L002_R2_001.blastx.daa](https://gannet.fish.washington.edu/Atumefaciens/20200103_cbai_diamond_blastx/329777_S4_L002_R2_001.blastx.daa) (31GB)
