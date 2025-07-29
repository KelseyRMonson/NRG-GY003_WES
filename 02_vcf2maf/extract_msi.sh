#!/bin/bash

output_file="msisensor_summary.txt"
base_path="/sc/arion/projects/NGSCRC/Work/Seqera/NRG-GY003_All_Varcalls/variant_calling/msisensorpro"

# Create header
echo -e "Sample_ID\tTumor_ID\tNormal_ID\tPercent" > "$output_file"

# Loop through each sample directory
for sample_dir in "$base_path"/*; do
    if [ -d "$sample_dir" ]; then
        sample_id=$(basename "$sample_dir")
        
        # Extract tumor and normal IDs from the sample name using _vs_ as delimiter
        tumor_id="${sample_id%%_vs_*}"
        normal_id="${sample_id#*_vs_}"
        
        # Look for the file with the same name as the directory (no extension)
        msi_file="$sample_dir/$sample_id"
        
        if [ -f "$msi_file" ]; then
            # Extract the third column (%) - skip header, create one row per data line
            awk -v sample="$sample_id" -v tumor="$tumor_id" -v normal="$normal_id" 'NR>1 {print sample "\t" tumor "\t" normal "\t" $3}' "$msi_file" >> "$output_file"
        else
            echo "Warning: MSI file not found for $sample_id" >&2
        fi
    fi
done

echo "Output written to $output_file"
