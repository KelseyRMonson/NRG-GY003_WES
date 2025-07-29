#!/bin/bash

module load bcftools

# Set output directory and log file
output_base_dir="/sc/arion/projects/NGSCRC/Work/NRG-GY003/WES/Consensus_VCFs"
output_log="${output_base_dir}/consensus_vcf_processing_log.txt"

# Create base output directory
mkdir -p "$output_base_dir"

echo "Consensus VCF Processing Log - $(date)" > "$output_log"
echo "=======================================" >> "$output_log"
echo "" >> "$output_log"

# Function to log messages to both console and file
log_message() {
    echo "$1"
    echo "$1" >> "$output_log"
}

# Function to check if sample is tumor-related
is_tumor_sample() {
    local sample_name="$1"
    
    # Check for tumor-only pattern (ends with -Tu-XXXXXX)
    if [[ "$sample_name" =~ -Tu-[A-Z0-9]+$ ]]; then
        return 0  # true - tumor-only sample
    fi
    
    # Check for tumor vs normal pattern (contains -Tu-XXXXXX_vs_)
    if [[ "$sample_name" =~ -Tu-[A-Z0-9]+_vs_ ]]; then
        return 0  # true - tumor vs normal sample
    fi
    
    return 1  # false - not a tumor sample
}

# Function to extract patient ID from sample name
get_patient_id() {
    local sample_name="$1"
    
    # For tumor-only samples like "GADCNC-Tu-0CDQZE", extract "GADCNC"
    if [[ "$sample_name" =~ ^([A-Z0-9]+)-Tu-[A-Z0-9]+$ ]]; then
        echo "${BASH_REMATCH[1]}"
        return
    fi
    
    # For tumor vs normal samples like "GABJHS-Tu-0CDQO6_vs_GABJHS-N-0C9S7D", extract "GABJHS"
    if [[ "$sample_name" =~ ^([A-Z0-9]+)-Tu-[A-Z0-9]+_vs_ ]]; then
        echo "${BASH_REMATCH[1]}"
        return
    fi
    
    # If no pattern matches, return the original name
    echo "$sample_name"
}

# Function to check if variant caller is for SNP/indels (not SVs)
is_snp_indel_caller() {
    local variantcaller="$1"
    local vcf_path="$2"
    
    # Exclude structural variant callers
    case "$variantcaller" in
        "tiddit"|"manta"|"delly"|"lumpy"|"gridss"|"svaba")
            return 1  # false - this is an SV caller
            ;;
        "mutect2"|"freebayes"|"strelka"|"haplotypecaller"|"varscan"|"bcftools")
            return 0  # true - these are SNP/indel callers
            ;;
        *)
            # Also check the filename for SV indicators
            if [[ "$vcf_path" =~ (sv|SV|structural) ]]; then
                return 1  # false - filename suggests SV
            fi
            return 0  # true - assume SNP/indel caller if not explicitly SV
            ;;
    esac
}

# Function to merge Strelka SNV and indel files if both exist
merge_strelka_files() {
    local sample="$1"
    local sample_dir="$2"
    local -n vcf_array=$3
    
    local strelka_snv=""
    local strelka_indel=""
    local other_vcfs=()
    
    # Separate Strelka files from others
    for vcf in "${vcf_array[@]}"; do
        if [[ "$vcf" =~ strelka.*somatic_snvs ]]; then
            strelka_snv="$vcf"
        elif [[ "$vcf" =~ strelka.*somatic_indels ]]; then
            strelka_indel="$vcf"
        else
            other_vcfs+=("$vcf")
        fi
    done
    
    # If we have both Strelka SNV and indel files, merge them
    if [ -n "$strelka_snv" ] && [ -n "$strelka_indel" ]; then
        log_message "Found Strelka SNV and indel files - merging them"
        
        local merged_strelka="${sample_dir}/strelka_merged_${sample}.vcf.gz"
        
        # Check if both files exist and are not empty
        if [ -f "$strelka_snv" ] && [ -s "$strelka_snv" ] && [ -f "$strelka_indel" ] && [ -s "$strelka_indel" ]; then
            # Create indexes if they don't exist
            if [ ! -f "${strelka_snv}.tbi" ]; then
                bcftools index -t "$strelka_snv"
            fi
            if [ ! -f "${strelka_indel}.tbi" ]; then
                bcftools index -t "$strelka_indel"
            fi
            
            # Sort both files first, then merge
            if bcftools sort "$strelka_snv" -O z -o "${sample_dir}/temp_snv_sorted.vcf.gz" && \
               bcftools sort "$strelka_indel" -O z -o "${sample_dir}/temp_indel_sorted.vcf.gz" && \
               bcftools index -t "${sample_dir}/temp_snv_sorted.vcf.gz" && \
               bcftools index -t "${sample_dir}/temp_indel_sorted.vcf.gz" && \
               bcftools concat "${sample_dir}/temp_snv_sorted.vcf.gz" "${sample_dir}/temp_indel_sorted.vcf.gz" -a -O z -o "$merged_strelka"; then
                
                bcftools index -t "$merged_strelka"
                log_message "Successfully merged Strelka files"
                
                # Clean up temporary files
                rm -f "${sample_dir}/temp_snv_sorted.vcf.gz" "${sample_dir}/temp_snv_sorted.vcf.gz.tbi" "${sample_dir}/temp_indel_sorted.vcf.gz" "${sample_dir}/temp_indel_sorted.vcf.gz.tbi"
                
                # Update the array to include merged file instead of separate files
                vcf_array=("${other_vcfs[@]}" "$merged_strelka")
            else
                log_message "ERROR: Failed to merge Strelka files, using separate files"
                # Clean up any partial temporary files
                rm -f "${sample_dir}/temp_snv_sorted.vcf.gz" "${sample_dir}/temp_snv_sorted.vcf.gz.tbi" "${sample_dir}/temp_indel_sorted.vcf.gz" "${sample_dir}/temp_indel_sorted.vcf.gz.tbi"
                vcf_array=("${other_vcfs[@]}" "$strelka_snv" "$strelka_indel")
            fi
        else
            log_message "WARNING: One or both Strelka files are missing/empty, using available files"
            if [ -f "$strelka_snv" ] && [ -s "$strelka_snv" ]; then
                other_vcfs+=("$strelka_snv")
            fi
            if [ -f "$strelka_indel" ] && [ -s "$strelka_indel" ]; then
                other_vcfs+=("$strelka_indel")
            fi
            vcf_array=("${other_vcfs[@]}")
        fi
    fi
}

# Function to generate consensus for a specific threshold (returns count only)
generate_consensus() {
    local sample="$1"
    local sample_dir="$2"
    local threshold="$3"
    local sites_file="$4"
    local indexed_vcfs=("${@:5}")
    
    local suffix=""
    if [ "$threshold" -eq 2 ]; then
        suffix="2plus_callers"
    elif [ "$threshold" -eq 3 ]; then
        suffix="3plus_callers"
    else
        suffix="${threshold}plus_callers"
    fi
    
    # Count variants that appear in at least the threshold number of files
    local consensus_count=$(awk -v thresh=$threshold '{
        caller_count = 0
        for(i=1; i<=length($5); i++) {
            if(substr($5,i,1) == "1") caller_count++
        }
        if(caller_count >= thresh) print
    }' "$sites_file" | wc -l)
    
    if [ "$consensus_count" -gt 0 ]; then
        # Extract positions that appear in at least threshold files
        local consensus_positions="${sample_dir}/consensus_positions_${suffix}.txt"
        awk -v thresh=$threshold '{
            caller_count = 0
            for(i=1; i<=length($5); i++) {
                if(substr($5,i,1) == "1") caller_count++
            }
            if(caller_count >= thresh) print $1"\t"$2
        }' "$sites_file" > "$consensus_positions"
        
        # Create consensus VCF using the first input file as template
        local consensus_vcf="${sample_dir}/consensus_${sample}_${suffix}.vcf.gz"
        if bcftools view -R "$consensus_positions" "${indexed_vcfs[0]}" -O z -o "$consensus_vcf" 2>/dev/null; then
            bcftools index -t "$consensus_vcf" 2>/dev/null
            
            # Count final variants in consensus VCF
            local final_variant_count=$(bcftools view -H "$consensus_vcf" | wc -l)
            
            # Log success message (redirect to avoid capture)
            {
                log_message "SUCCESS: Consensus VCF (=${threshold} callers) created: $(basename "$consensus_vcf")"
                log_message "  Final consensus variants: $final_variant_count"
                log_message "Consensus VCF statistics (=${threshold} callers):"
                bcftools stats "$consensus_vcf" | grep "^SN" >> "$output_log"
            } >&2
        else
            {
                log_message "ERROR: Failed to create consensus VCF for =${threshold} callers"
            } >&2
            consensus_count=0
        fi
        
        rm "$consensus_positions" 2>/dev/null
    else
        {
            log_message "RESULT: No consensus variants found for =${threshold} callers"
        } >&2
    fi
    
    # Return only the consensus count (no other output)
    echo "$consensus_count"
}

# Function to process a single sample
process_sample() {
    sample=$1
    shift
    vcf_files=("$@")
    
    # Create sample-specific directory
    sample_dir="${output_base_dir}/${sample}"
    mkdir -p "$sample_dir"
    
    log_message "Processing tumor sample: $sample"
    log_message "Output directory: $sample_dir"
    
    # Merge Strelka files if needed
    merge_strelka_files "$sample" "$sample_dir" vcf_files
    
    log_message "Number of VCF files: ${#vcf_files[@]}"
    
    # Count total variants in input files
    total_input_variants=0
    for vcf in "${vcf_files[@]}"; do
        if [ -f "$vcf" ] && [ -s "$vcf" ]; then
            variant_count=$(bcftools view -H "$vcf" | wc -l)
            log_message "  $(basename "$vcf"): $variant_count variants"
            total_input_variants=$((total_input_variants + variant_count))
        fi
    done
    log_message "Total input variants across all SNP/indel callers: $total_input_variants"
    
    # Check if all VCF files exist and create indexes if needed
    indexed_vcfs=()
    for vcf in "${vcf_files[@]}"; do
        if [ ! -f "$vcf" ]; then
            log_message "ERROR: VCF file $vcf does not exist"
            return 1
        fi
        if [ ! -s "$vcf" ]; then
            log_message "WARNING: VCF file $vcf is empty"
            continue
        fi
        
        # Check if index exists, create if not
        index_file="${vcf}.tbi"
        if [ ! -f "$index_file" ]; then
            if ! bcftools index -t "$vcf"; then
                log_message "ERROR: Failed to create index for $vcf"
                return 1
            fi
        fi
        
        indexed_vcfs+=("$vcf")
    done
    
    if [ "${#indexed_vcfs[@]}" -lt 2 ]; then
        log_message "ERROR: Need at least 2 valid SNP/indel VCF files for consensus, found ${#indexed_vcfs[@]}"
        return 1
    fi
    
    log_message "Running consensus analysis for ${#indexed_vcfs[@]} SNP/indel callers"
    
    # Run bcftools isec to get all intersection statistics
    isec_dir="${sample_dir}/isec_all"
    mkdir -p "$isec_dir"
    
    if ! bcftools isec -O v -p "$isec_dir" "${indexed_vcfs[@]}"; then
        log_message "ERROR: bcftools isec failed for sample $sample"
        rm -rf "$isec_dir"
        return 1
    fi
    
    # Analyze the complete intersection results
    sites_file="${isec_dir}/sites.txt"
    if [ -f "$sites_file" ]; then
        # Count total unique positions
        total_unique_positions=$(wc -l < "$sites_file")
        
        # Count variants by caller number using the binary string
        declare -A caller_counts
        for i in $(seq 1 ${#indexed_vcfs[@]}); do
            caller_counts[$i]=$(awk -v count=$i '{
                caller_count = 0
                for(j=1; j<=length($5); j++) {
                    if(substr($5,j,1) == "1") caller_count++
                }
                if(caller_count == count) print
            }' "$sites_file" | wc -l)
        done
        
        log_message "Variant filtering summary for $sample:"
        log_message "  Total unique variant positions: $total_unique_positions"
        log_message "  Variants by caller count:"
        for i in $(seq 1 ${#indexed_vcfs[@]}); do
            log_message "    - Found in $i caller(s): ${caller_counts[$i]}"
        done
        
        # Generate consensus files for different thresholds
        if [ "${#indexed_vcfs[@]}" -ge 2 ]; then
            # Generate 2+ caller consensus
            consensus_2plus=$(generate_consensus "$sample" "$sample_dir" 2 "$sites_file" "${indexed_vcfs[@]}" 2>/dev/null)
            filtered_out_2plus=$(awk '{
                caller_count = 0
                for(i=1; i<=length($5); i++) {
                    if(substr($5,i,1) == "1") caller_count++
                }
                if(caller_count < 2) print
            }' "$sites_file" | wc -l)
            
            log_message "  Consensus variants (=2 callers): $consensus_2plus"
            log_message "  Total filtered out (=2 callers): $filtered_out_2plus"
            if [ "$total_unique_positions" -gt 0 ]; then
                retention_rate_2plus=$(awk "BEGIN {printf \"%.2f\", $consensus_2plus/$total_unique_positions*100}")
                log_message "  Retention rate (=2 callers): ${retention_rate_2plus}%"
            fi
        fi
        
        if [ "${#indexed_vcfs[@]}" -ge 3 ]; then
            # Generate 3+ caller consensus
            consensus_3plus=$(generate_consensus "$sample" "$sample_dir" 3 "$sites_file" "${indexed_vcfs[@]}" 2>/dev/null)
            filtered_out_3plus=$(awk '{
                caller_count = 0
                for(i=1; i<=length($5); i++) {
                    if(substr($5,i,1) == "1") caller_count++
                }
                if(caller_count < 3) print
            }' "$sites_file" | wc -l)
            
            log_message "  Consensus variants (=3 callers): $consensus_3plus"
            log_message "  Total filtered out (=3 callers): $filtered_out_3plus"
            if [ "$total_unique_positions" -gt 0 ]; then
                retention_rate_3plus=$(awk "BEGIN {printf \"%.2f\", $consensus_3plus/$total_unique_positions*100}")
                log_message "  Retention rate (=3 callers): ${retention_rate_3plus}%"
            fi
        fi
        
    else
        log_message "ERROR: No sites.txt file found - bcftools isec may have failed"
    fi
    
    # Clean up temporary merged files and isec directory
    rm -f "${sample_dir}/strelka_merged_${sample}.vcf.gz" "${sample_dir}/strelka_merged_${sample}.vcf.gz.tbi"
    rm -rf "$isec_dir"
    
    log_message "Completed processing for sample: $sample"
    log_message "---"
}

# Main script
csv_dir="/sc/arion/projects/NGSCRC/Work/Seqera/NRG-GY003_All_Varcalls/csv"
input_csv="$csv_dir/variantcalled.csv"

log_message "Starting consensus VCF processing for TUMOR samples only (SNP/indel callers including Mutect2)"
log_message "Input CSV: $input_csv"
log_message "Output directory: $output_base_dir"

# Check if input file exists
if [ ! -f "$input_csv" ]; then
    log_message "Error: Input file $input_csv not found!"
    exit 1
fi

# Initialize variables
declare -A sample_vcfs  # Associative array to group VCFs by patient ID
declare -A sample_names # Map patient ID to full sample name
total_samples=0
successful_samples=0
skipped_normal_samples=0
skipped_insufficient_vcfs=0
skipped_sv_calls=0

# Skip header and process each sample
while IFS=',' read -r patient sample variantcaller vcf; do
    # Remove any quotes and whitespace
    patient=$(echo "$patient" | tr -d '"' | xargs)
    sample=$(echo "$sample" | tr -d '"' | xargs)
    variantcaller=$(echo "$variantcaller" | tr -d '"' | xargs)
    vcf=$(echo "$vcf" | tr -d '"' | xargs)
    
    # Check if VCF file exists
    if [ ! -f "$vcf" ]; then
        continue
    fi
    
    # Skip structural variant callers
    if ! is_snp_indel_caller "$variantcaller" "$vcf"; then
        skipped_sv_calls=$((skipped_sv_calls + 1))
        continue
    fi
    
    # Determine the key to group by
    if is_tumor_sample "$sample"; then
        # For tumor samples, extract patient ID and use that as the key
        patient_id=$(get_patient_id "$sample")
        sample_key="$patient_id"
        sample_names["$patient_id"]="$sample"
    elif [[ "$variantcaller" == "mutect2" ]]; then
        # For Mutect2, the sample name is likely the patient ID
        sample_key="$sample"
        sample_names["$sample"]="$sample (Mutect2)"
    else
        # Skip non-tumor, non-Mutect2 samples
        skipped_normal_samples=$((skipped_normal_samples + 1))
        continue
    fi
    
    # Add VCF to the appropriate group
    if [[ -n "${sample_vcfs[$sample_key]}" ]]; then
        sample_vcfs["$sample_key"]+=" $vcf"
    else
        sample_vcfs["$sample_key"]="$vcf"
    fi
    
done < <(tail -n +2 "$input_csv")

# Process each patient's VCFs
for patient_id in "${!sample_vcfs[@]}"; do
    # Convert space-separated string to array
    IFS=' ' read -ra vcf_files <<< "${sample_vcfs[$patient_id]}"
    sample_name="${sample_names[$patient_id]}"
    
    if [ "${#vcf_files[@]}" -ge 2 ]; then
        total_samples=$((total_samples + 1))
        if process_sample "$patient_id" "${vcf_files[@]}"; then
            successful_samples=$((successful_samples + 1))
        fi
    else
        log_message "Skipping $patient_id: fewer than 2 SNP/indel VCFs (found ${#vcf_files[@]}) - need =2 for consensus"
        skipped_insufficient_vcfs=$((skipped_insufficient_vcfs + 1))
    fi
done

# Final summary
log_message ""
log_message "======================================="
log_message "FINAL SUMMARY"
log_message "======================================="
log_message "Total tumor samples processed: $total_samples"
log_message "Successful tumor samples: $successful_samples"
log_message "Failed tumor samples: $((total_samples - successful_samples))"
log_message "Skipped normal samples: $skipped_normal_samples"
log_message "Skipped tumor samples (insufficient SNP/indel VCFs): $skipped_insufficient_vcfs"
log_message "Skipped SV variant calls: $skipped_sv_calls"
log_message ""
log_message "Output files created:"
for sample_dir in "${output_base_dir}"/*/; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        log_message "Sample: $sample_name"
        for consensus_file in "$sample_dir"/consensus_*.vcf.gz; do
            if [ -f "$consensus_file" ]; then
                file_size=$(ls -lh "$consensus_file" | awk '{print $5}')
                variant_count=$(bcftools view -H "$consensus_file" | wc -l)
                log_message "  $(basename "$consensus_file") (Size: $file_size, Variants: $variant_count)"
            fi
        done
    fi
done

log_message ""
log_message "Processing completed at: $(date)"
log_message "Log file: $output_log"

echo ""
echo "Processing complete! Check the log file: $output_log"
echo "Output directory: $output_base_dir"
echo "Total samples processed: $successful_samples"