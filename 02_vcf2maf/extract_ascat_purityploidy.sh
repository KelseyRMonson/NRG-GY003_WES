#!/bin/bash

output_file="ascat_summary.txt"
base_path="/sc/arion/projects/NGSCRC/Work/Seqera/NRG-GY003_All_Varcalls/variant_calling/ascat"

# Create header
echo -e "Sample_ID\tTumor_ID\tNormal_ID\tPurity\tPloidy" > "$output_file"

# Loop through each sample directory
for sample_dir in "$base_path"/*; do
    if [ -d "$sample_dir" ]; then
        sample_id=$(basename "$sample_dir")
        
        # Extract tumor and normal IDs from the sample name using _vs_ as delimiter
        tumor_id="${sample_id%%_vs_*}"
        normal_id="${sample_id#*_vs_}"
        
        purity_file="$sample_dir/${sample_id}.purityploidy.txt"
        
        if [ -f "$purity_file" ]; then
            # Skip header, extract AberrantCellFraction (purity) and Ploidy from the second line
            purity=$(awk 'NR==2 {print $1}' "$purity_file")
            ploidy=$(awk 'NR==2 {print $2}' "$purity_file")
            echo -e "$sample_id\t$tumor_id\t$normal_id\t$purity\t$ploidy" >> "$output_file"
        else
            echo -e "$sample_id\t$tumor_id\t$normal_id\tNA\tNA" >> "$output_file"
            echo "Warning: Purity file not found for $sample_id" >&2
        fi
    fi
done

echo "Output written to $output_file"
