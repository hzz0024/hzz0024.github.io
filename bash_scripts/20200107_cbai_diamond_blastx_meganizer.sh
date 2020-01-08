#!/bin/bash

# Script to run MEGAN6 meganizer on DIAMOND DAA files from
# 20200103_cbai_diamond_blastx Mox job.

# Requires MEGAN mapping files from:
# http://ab.inf.uni-tuebingen.de/data/software/megan6/download/

# Program path
meganizer=/home/sam/programs/megan/tools/daa-meganizer

# MEGAN mapping files
prot_acc2tax=/home/sam/data/databases/MEGAN/prot_acc2tax-Jul2019X1.abin
acc2interpro=/home/sam/data/databases/MEGAN/acc2interpro-Jul2019X.abin
acc2eggnog=/home/sam/data/databases/MEGAN/acc2eggnog-Jul2019X.abin

# Variables
threads=20

## Run MEGANIZER

# Capture start "time"
start=${SECONDS}
for daa in *.daa
do
  ${meganizer} \
  --in "${daa}" \
	--threads "${threads}" \
	--acc2taxa ${prot_acc2tax} \
	--acc2interpro2go ${acc2interpro} \
	--acc2eggnog ${acc2eggnog}
done

# Caputure end "time"
end=${SECONDS}

runtime=$((end-start))

# Print MEGANIZER runtime, in seconds
echo "Runtime was: ${runtime} seconds"