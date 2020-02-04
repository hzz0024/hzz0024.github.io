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
#SBATCH --chdir=/gscratch/scrubbed/samwhite/outputs/20200128_cbai_DEG_inf-vs-uninf

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

prefix="${timestamp}.${species}"
fasta_prefix="20200122.C_bairdi.megan.Trinity"

# Declare empty arrays
fastq_R1_array=()
fastq_R2_array=()

declare -A read_pairs_array

# Create comparisons array
comparisons=(
"infected-vs-uninfected" \
"D12-vs-D26" \
"D12_infected-vs-D12_uninfected" \
"D12_uninfected-vs-D26_uninfected" \
"D12_infected-vs-D26_infected" \
"D26_infected-vs-D26_uninfected"
)

# Create associative arrays

## Infection status
declare -A inf_status_array=( [329774]=infected [329775]=uninfected [329776]=infected [329777]=uninfected )
## Sampling day
declare -A sample_day_array=( [329774]=D12 [329775]=D12 [329776]=D26 [329777]=D26 )

for fastq in ${trimmed_reads_dir}/20200131*R1*.fq
do
	fastq_R1_array+=(${fastq})
done

# Functions
## Sample list file formatting
sample_list () {
	printf "%s\t%s\t%s\t%s\n" "${inf_status}" "${inf_status}_${sample_day}_0${counter}" "${fastq}" "${read_pairs_array[fastq]}" \
	>> "${comparison}.samples.txt"
}


for fastq in ${trimmed_reads_dir}/20200131*R2*.fq
do
	fastq_R2_array+=(${fastq})
done

for index in "${!fastq_R1_array[@]}"
do
	R1=${fastq_R1_array[index]}
	R2=${fastq_R2_array[index]}
	read_pairs_array+=([$R1]=$R2)
done


# Comparisons

for comparison in ${!comparisons[@]}
do
	counter=0

	# Record number of fiels in comparison
	field_count=$(echo ${comparisons[${comparison}]} | awk -F[_-] '{print NF}')

	# Check day/infection status
	# Comparison ayout is always day, inf, vs, day, inf
	inf_check=$(echo ${comparisons[${comparison}]} | awk -F[_-] '{print $1}')
	day_check1=$(echo ${comparisons[${comparison}]} | awk -F[_-] '{print $1}')
	day_check2=$(echo ${comparisons[${comparison}]} | awk -F[_-] '{print $4}')
	inf_check1=$(echo ${comparisons[${comparison}]} | awk -F[_-] '{print $2}')
	inf_check2=$(echo ${comparisons[${comparison}]} | awk -F[_-] '{print $5}')

	for fastq in ${!read_pairs_array[@]}
  do
  	fastq_nopath=${fastq##*/}
  	sample=$(echo ${fastq_nopath} | awk -F "." '{print $3}')
		inf_status=${inf_status_array[${sample}]}
		sample_day=${sample_day_array[$sample]}
    # If the field_count is equal to 3
  	# then we know the comparison is either inf_vs_uninf or d12_vs_d26
  	if [[ ${field_count} -eq 3 ]]; then
  		# If the inf_check string matches "infected", then print the associated values to sample file
  	  if [[ "${inf_check}" == "infected" ]]; then
  			if [[ "${sample}" = "329774" ]] || [[ "${sample}" = "329775" ]] || [[ "${sample}" = "329776" ]] || [[ "${sample}" = "329777" ]]; then
  				(( counter ++ ))
  	  		sample_list
  			else
  				if [[ "${sample}" = "329774" ]] || [[ "${sample}" = "329775" ]] || [[ "${sample}" = "329776" ]] || [[ "${sample}" = "329777" ]]; then
  					(( counter ++ ))
  		  		sample_list
  				fi
  			fi
  	  fi
  	else
  		if [[ "${day_check1}" == "D12" ]] && [[ "${day_check2}" == "D12" ]]; then
  			#statements
  			if [[ "${sample}" = "329774" ]] || [[ "${sample}" = "329775" ]] || [[ "${sample}" = "329776" ]] || [[ "${sample}" = "329777" ]]; then
  	  		(( counter ++ ))
  	  		sample_list
  			elif [[ "${day_check1}" == "D12" ]] && [[ "${day_check2}" == "D26" ]] && [[ "${inf_check1}" == "uninfected" ]] && [[ "${inf_check2}" == "uninfected" ]]; then
  				#statements
  				if [[ "${sample}" = "329775" ]] || [[ "${sample}" = "329777" ]]; then
  		  		(( counter ++ ))
  		  		sample_list
  				fi
  			elif [[ "${day_check1}" == "D12" ]] && [[ "${day_check2}" == "D26" ]] && [[ "${inf_check1}" == "infected" ]] && [[ "${inf_check2}" == "infected" ]]; then
  				#statements
  				if [[ "${sample}" = "329774" ]] || [[ "${sample}" = "329776" ]]; then
  		  		(( counter ++ ))
  		  		sample_list
  				fi
  			else
  				if [[ "${sample}" = "329776" ]] || [[ "${sample}" = "329777" ]]; then
  			  		(( counter ++ ))
  			  		sample_list
  				fi
  		fi
  	fi
  done
done



## Set input file locations
trimmed_reads_dir="/gscratch/srlab/sam/data/C_bairdi/RNAseq"
salmon_out_dir=""
transcriptome_dir="/gscratch/srlab/sam/data/C_bairdi/transcriptomes"
transcriptome="${transcriptome_dir}/${fasta_prefix}.fasta"
fasta_index="${transcriptome_dir}/${fasta_prefix}.fasta.fai"
fasta_seq_lengths="${transcriptome_dir}/${fasta_prefix}.seq_lens"
samples=""

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
mv ${trimmed_reads_dir}/[BN]* \
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
