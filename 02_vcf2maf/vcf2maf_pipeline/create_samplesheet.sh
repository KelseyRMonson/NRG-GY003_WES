#!/bin/bash

# Read your original samplelist.txt and create CSV
echo "sample_id,somatic_vcf,germline_vcf,has_germline" > samplesheet.csv

while IFS= read -r SAMPLE_PATH || [ -n "$SAMPLE_PATH" ]; do
    # Skip empty lines and comments
    if [[ -z "$SAMPLE_PATH" || "$SAMPLE_PATH" =~ ^[[:space:]]*# ]]; then
        continue
    fi
    
    # Remove leading/trailing whitespace
    SAMPLE_PATH=$(echo "$SAMPLE_PATH" | xargs)
    SAMPLE_ID=$(basename "$SAMPLE_PATH")
    
    SOMATIC_VCF="${SAMPLE_PATH}/final_annotated_refined_high_quality_consensus_${SAMPLE_ID}_2plus_callers_SOMATIC.vcf.gz"
    GERMLINE_VCF="${SAMPLE_PATH}/final_consensus_pathogenic_germline_${SAMPLE_ID}.vcf.gz"
    
    # Check if files exist
    if [[ -f "$SOMATIC_VCF" ]]; then
        if [[ -f "$GERMLINE_VCF" ]]; then
            # Both somatic and germline exist
            echo "${SAMPLE_ID},${SOMATIC_VCF},${GERMLINE_VCF},true" >> samplesheet.csv
        else
            # Only somatic exists (tumor-only)
            echo "${SAMPLE_ID},${SOMATIC_VCF},,false" >> samplesheet.csv
        fi
    else
        echo "Warning: Somatic VCF not found for $SAMPLE_ID: $SOMATIC_VCF"
    fi
done < samplelist.txt

# Print summary
TOTAL_SAMPLES=$(tail -n +2 samplesheet.csv | wc -l)
PAIRED_SAMPLES=$(tail -n +2 samplesheet.csv | awk -F',' '$4=="true"' | wc -l)
TUMOR_ONLY_SAMPLES=$(tail -n +2 samplesheet.csv | awk -F',' '$4=="false"' | wc -l)

echo "Created samplesheet.csv with:"
echo "  Total samples: $TOTAL_SAMPLES"
echo "  Paired samples (tumor + germline): $PAIRED_SAMPLES"
echo "  Tumor-only samples: $TUMOR_ONLY_SAMPLES"