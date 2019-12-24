#!/bin/bash
## Job Name
#SBATCH --job-name=blastx_cbai
## Allocation Definition
#SBATCH --account=coenv
#SBATCH --partition=coenv
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=25-00:00:00
## Memory per node
#SBATCH --mem=120G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samwhite@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/samwhite/outputs/20191224_cbai_blastx_outfmt-11

# Load Python Mox module for Python module availability

module load intel-python3_2017

# Document programs in PATH (primarily for program version ID)

{
date
echo ""
echo "System PATH for $SLURM_JOB_ID"
echo ""
printf "%0.s-" {1..10}
echo "${PATH}" | tr : \\n
} >> system_path.log


wd="$(pwd)"
timestamp=$(date +%Y%m%d)

# Paths to input/output files
blastx_out="${wd}/${timestamp}-20191218.C_bairdi.Trinity.fasta.blastx.asn"
sp_db="/gscratch/srlab/programs/Trinotate-v3.1.1/admin/uniprot_sprot.pep"

trinity_fasta="/gscratch/scrubbed/samwhite/outputs/20191218_cbai_trinity_RNAseq/trinity_out_dir/20191218.C_bairdi.Trinity.fasta"

# Paths to programs
blast_dir="/gscratch/srlab/programs/ncbi-blast-2.8.1+/bin"
blastx="${blast_dir}/blastx"

threads=28

# Run blastx on Trinity fasta
"${blastx}" \
-query "${trinity_fasta}" \
-db "${sp_db}" \
-max_target_seqs 1 \
-outfmt 11 \
-evalue 1e-4 \
-num_threads "${threads}" \
> "${blastx_out}"
