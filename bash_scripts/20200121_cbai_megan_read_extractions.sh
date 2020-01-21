



# Paths to data
crab_data=/home/sam/data/C_bairdi/RNAseq
hemat_data=/home/sam/data/Hematodinium/RNAseq

# Paths to programs
samtools=/path/to/samtools
seqtk=/path/to/seqtk


timestamp=$(date +%Y%m%d)


for directory in ${crab_data} ${hemat_data}
do
	# Get species name
	species=$(echo ${directory} | awk -F"/" '{print $5}')
	cd ${directory}

	######################################################
	# Create FastA IDs list to use for sequence extraction
	######################################################
	# Index FastA files
	for fasta in *.fasta
	do
		${samtools} faidx ${fasta}
	done

  # Extract sequence IDs
	for fai in *.fai
	do
		awk '{print $1}' >> ${timestamp}.${species}.seqtk.read_id.list
	done
done

# Check to ensure that read_id.list files only have unique FastA IDs

# Compare FastA and FastA index counts
wc -l 20200121_cbai_304428-reads-Alveolata.fasta.fai
grep ">" -c 20200121_cbai_304428-reads-Alveolata.fasta | wc -l

# Check that all FastA IDs are unique
awk '{print $1}' 20200121_cbai_304428-reads-Alveolata.fasta.fai | head -n 5
awk '{print $1}' 20200121_cbai_304428-reads-Alveolata.fasta.fai | sort | uniq | wc -l




# Create array of fastq R1 files
for fastq in *R1*.gz
do
  ${seqtk} subseq ${fastq} seqtk.id.list >> ${timestamp}.megan_R1.fq
  fastq_array_R1+=("${fastq}")
done

# Create array of fastq R2 files
for fastq in *R2*.gz
do
  ${seqtk} subseq ${fastq} seqtk.id.list >> ${timestamp}.megan_R2.fq
  fastq_array_R2+=("${fastq}")
done
