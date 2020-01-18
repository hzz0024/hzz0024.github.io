#!/bin/bash

# Script to run MEGAN6 meganizer on DIAMOND DAA files from
# 20200103_cbai_diamond_blastx Mox job.

# Requires MEGAN mapping files from:
# http://ab.inf.uni-tuebingen.de/data/software/megan6/download/

# Exit script if any command fails
set -e

# Program path
meganizer=/home/sam/programs/megan/tools/daa2rma

# MEGAN mapping files
prot_acc2tax=/home/sam/data/databases/MEGAN/prot_acc2tax-Jul2019X1.abin
acc2interpro=/home/sam/data/databases/MEGAN/acc2interpro-Jul2019X.abin
acc2eggnog=/home/sam/data/databases/MEGAN/acc2eggnog-Jul2019X.abin


## Inititalize arrays
l001_array_R1=()
l001_array_R2=()
l002_array_R1=()
l002_array_R2=()
daa_array_R1=()
daa_array_R2=()

# Create arrays for sequencing lanes
for daa in *L001_R1*
do
  l001_array_R1+=("${daa}")
done

for daa in *L001_R2*
do
  l001_array_R2+=("${daa}")
done

for daa in *l002_array_R1*
do
  l002_array_R1+=("${daa}")
done

for daa in *l002_array_R2*
do
  l002_array_R2+=("${daa}")
done



# Concatenate R1 and R2 DAA files
for index in "${!l001_array_R1[@]}"
do
  sample_name=$(echo "${l001_array_R1[index]}" | awk -F "_" '{print $1}')
  #  Concatenate R1 DAA files
  {
    cat "${l001_array_R1[index]}"
    cat "${l002_array_R1[index]}"
  } >> "${sample_name}"R1_.blastx.cat.daa

  # Concatenate R2 DAA files
  {
    cat "${l001_array_R2[index]}"
    cat "${l002_array_R2[index]}"
  } >> "${sample_name}"R2_.blastx.cat.daa
done

# Create array of DAA R1 files
for daa in *R1*.cat.daa
do
  daa_array_R1+=("${daa}")
done

# Create array of DAA R2 files
for daa in *R2*.cat.daa
do
  daa_array_R2+=("${daa}")
done

## Run MEGANIZER

# Capture start "time"
start=${SECONDS}
for index in "${!daa_array_R1[@]}"
do
  sample_name=$(echo "${daa_array_R1[index]}" | awk -F "_" '{print $1}')

  # Run daa2rma with paired option
  ${meganizer} \
  --paired \
  --in "${daa_array_R1[index]}" "${daa_array_R2[index]}" \
	--acc2taxa ${prot_acc2tax} \
	--acc2interpro2go ${acc2interpro} \
	--acc2eggnog ${acc2eggnog} \
  --out "${sample_name}".daa2rma.rma6 \
  2>&1 | tee --append daa2rma_log.txt
done

# Caputure end "time"
end=${SECONDS}

runtime=$((end-start))

# Print MEGANIZER runtime, in seconds

{
  echo ""
  echo "---------------------"
  echo ""
  echo "Total runtime was: ${runtime} seconds"
} >> daa2rma_log.txt
