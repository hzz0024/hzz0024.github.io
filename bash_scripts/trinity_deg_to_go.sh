
# Input file
## Expects Trinity edgeR GO enrichment format:
## category	over_represented_pvalue	under_represented_pvalue	numDEInCat	numInCat	term	ontology	over_represented_FDR	go_term	gene_ids
## Field 10 (gene_ids) contains comma separated gene_ids that fall in the given GO term in the "category" column
goseq="salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.infected-UP.subset.GOseq.enriched"

# Get GO term and gene_ids
awk -F "\t" '{print $1, $10}' \
${goseq} \
# Remove commas
| sed 's/,//g' \
# Create tab-delimited file by converting spaces to tabs
| tr ' ' $'\t' \
> $file

# Identify the first line number which contains a gene_id
begin_goterms=$(grep "TRINITY" "$file" | awk '{for (i=1;i<=NF;i++) if($i ~/TRINITY/) print i}' | sort -ug | head -n1)

# "Unfolds" gene_ids to a single gene_id per row
while read -r line
do
	# Capture the length of the longest row
	max_field=$(echo "$line" | awk -F "\t" '{print NF}')

	# Retain the GO term in the first field
	fixed_fields=$(echo "$line" | cut -f1)

	# Since not all the lines contain the same number of fields (e.g. may not have GO terms),
	# evaluate the number of fields in each line to determine how to handle current line.

	# If the value in max_field is less than the field number where the GO terms begin,
	# then just print the current line (%s) followed by a newline (\n).
	if (( "$max_field" < "$begin_goterms" ))
	then
		printf "%s\n" "$line"
	else goterms=$(echo "$line" | cut -f"$begin_goterms"-"$max_field")

  # Assign values in the variable "goterms" to a new indexed array (called "array"),
  # with tab delimiter (IFS=$'\t')
  IFS=$'\t' read -r -a array <<<"$goterms"

  # Iterate through each element of the array.
  # Print the first n fields (i.e. the fields stored in "fixed_fields") followed by a tab (%s\t).
  # Print the current element in the array (i.e. the current GO term) followed by a new line (%s\n).
  for element in "${!array[@]}"
  do
	  printf "%s\t%s\n" "$fixed_fields" "${array[$element]}"
  done
  fi
done < "$file" > "${output_file}"


while read line
do grep ${line} ${goseq}

done < goterms.txt
