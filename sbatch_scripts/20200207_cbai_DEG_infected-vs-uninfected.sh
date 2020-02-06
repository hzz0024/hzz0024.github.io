#!/bin/bash
## Job Name
#SBATCH --job-name=DEG_cbai_inf-vs-uninf
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=05-00:00:00
## Memory per node
#SBATCH --mem=120G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samwhite@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/samwhite/outputs/20200207_cbai_DEG/infected-vs-uninfected

# Exit script if any command fails
set -e

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
species="cbai"
threads=28
comparison="${wd##*/}"

prefix="${timestamp}.${species}"
fasta_prefix="20200122.C_bairdi.megan.Trinity"


## Set input file locations
trimmed_reads_dir="/gscratch/srlab/sam/data/C_bairdi/RNAseq"
salmon_out_dir="${wd}"
transcriptome_dir="/gscratch/srlab/sam/data/C_bairdi/transcriptomes"
transcriptome="${transcriptome_dir}/${fasta_prefix}.fasta"
fasta_index="${transcriptome_dir}/${fasta_prefix}.fasta.fai"
fasta_seq_lengths="${transcriptome_dir}/${fasta_prefix}.seq_lens"
samples="${comparison}.samples.txt"

gene_map="${transcriptome_dir}/${fasta_prefix}.gene_trans_map"
salmon_gene_matrix="${salmon_out_dir}/salmon.gene.TMM.EXPR.matrix"
salmon_iso_matrix="${salmon_out_dir}/salmon.isoform.TMM.EXPR.matrix"
go_annotations="${transcriptome_dir}/20200126.cbai.trinotate.go_annotations.txt"


# Standard output/error files
diff_expr_stdout="diff_expr_stdout.txt"
diff_expr_stderr="diff_expr_stderr.txt"
matrix_stdout="matrix_stdout.txt"
matrix_stderr="matrix_stderr.txt"
salmon_stdout="salmon_stdout.txt"
salmon_stderr="salmon_stderr.txt"
tpm_length_stdout="tpm_length_stdout.txt"
tpm_length_stderr="tpm_length_stderr.txt"
trinity_DE_stdout="trinity_DE_stdout.txt"
trinity_DE_stderr="trinity_DE_stderr.txt"

edgeR_dir=""

#programs
trinity_abundance=/gscratch/srlab/programs/trinityrnaseq-v2.9.0/util/align_and_estimate_abundance.pl
trinity_matrix=/gscratch/srlab/programs/trinityrnaseq-v2.9.0/util/abundance_estimates_to_matrix.pl
trinity_DE=/gscratch/srlab/programs/trinityrnaseq-v2.9.0/Analysis/DifferentialExpression/run_DE_analysis.pl
diff_expr=/gscratch/srlab/programs/trinityrnaseq-v2.9.0/Analysis/DifferentialExpression/analyze_diff_expr.pl
trinity_tpm_length=/gscratch/srlab/programs/trinityrnaseq-v2.9.0/util/misc/TPM_weighted_gene_length.py


cd ${trimmed_reads_dir}
time ${trinity_abundance} \
--output_dir ${salmon_out_dir} \
--transcripts ${transcriptome} \
--seqType fq \
--samples_file ${samples} \
--SS_lib_type RF \
--est_method salmon \
--aln_method bowtie2 \
--gene_trans_map "${gene_map}" \
--prep_reference \
--thread_count "${threads}" \
1> ${salmon_out_dir}/${salmon_stdout} \
2> ${salmon_out_dir}/${salmon_stderr}

# Move output folders
mv ${trimmed_reads_dir}/[iu]* \
${salmon_out_dir}

cd ${salmon_out_dir}

# Convert abundance estimates to matrix
${trinity_matrix} \
--est_method salmon \
--gene_trans_map ${gene_map} \
--out_prefix salmon \
--name_sample_by_basedir \
# NEED directcory/quant.sf - directory comes from samples list 2nd column
1> ${matrix_stdout} \
2> ${matrix_stderr}

# Generate weighted gene lengths
"${trinity_tpm_length}" \
--gene_trans_map "${gene_map}" \
--trans_lengths "${fasta_seq_lengths}" \
--TPM_matrix "${salmon_iso_matrix}" \
> Trinity.gene_lengths.txt \
2> ${tpm_length_stderr}

# Differential expression analysis
cd ${transcriptome_dir}
${trinity_DE} \
--matrix ${salmon_out_dir}/salmon.gene.counts.matrix \
--method edgeR \
--samples_file ${samples} \
1> ${trinity_DE_stdout} \
2> ${trinity_DE_stderr}

mv edgeR* ${salmon_out_dir}


# Run differential expression on edgeR output matrix
# Set fold difference to 2-fold (ie. -C 1 = 2^1)
# P value <= 0.05
# Has to run from edgeR output directory

# Pulls edgeR directory name and removes leading ./ in find output
cd ${salmon_out_dir}
edgeR_dir=$(find . -type d -name "edgeR*" | sed 's%./%%')
cd "${edgeR_dir}"
mv "${transcriptome_dir}/${trinity_DE_stdout}" .
mv "${transcriptome_dir}/${trinity_DE_stderr}" .
${diff_expr} \
--matrix "${salmon_gene_matrix}" \
--samples ${samples} \
--examine_GO_enrichment \
--GO_annots "${go_annotations}" \
--include_GOplot \
--gene_lengths ${salmon_out_dir}/Trinity.gene_lengths.txt \
-C 1 \
-P 0.05 \
1> ${diff_expr_stdout} \
2> ${diff_expr_stderr}
