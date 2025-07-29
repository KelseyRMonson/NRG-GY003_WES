#!/bin/bash

module load bcftools

# Set output directories
output_base_dir="/sc/arion/projects/NGSCRC/Work/NRG-GY003/WES/Consensus_VCFs"
high_quality_dir="/sc/arion/projects/NGSCRC/Work/NRG-GY003/WES/High_Quality_Consensus_VCFs"
output_log="${output_base_dir}/consensus_vcf_processing_log.txt"

# Create output directories
mkdir -p "$output_base_dir"
mkdir -p "$high_quality_dir"

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

# Function to detect caller type from VCF path or content
detect_caller_type() {
    local vcf_path="$1"
    local basename=$(basename "$vcf_path")
    
    # Check filename patterns first
    if [[ "$basename" =~ freebayes ]]; then
        echo "freebayes"
    elif [[ "$basename" =~ mutect2 ]]; then
        echo "mutect2"
    elif [[ "$basename" =~ strelka ]]; then
        echo "strelka"
    else
        # Check VCF header for caller information
        local header_info=$(bcftools view -h "$vcf_path" | grep -i "source\|command\|program" | head -5)
        if echo "$header_info" | grep -qi "freebayes"; then
            echo "freebayes"
        elif echo "$header_info" | grep -qi "mutect"; then
            echo "mutect2"
        elif echo "$header_info" | grep -qi "strelka"; then
            echo "strelka"
        else
            echo "unknown"
        fi
    fi
}

# Function to apply FreeBayes-specific quality filters
apply_freebayes_filters() {
    local input_vcf="$1"
    local output_vcf="$2"
    
    log_message "    Applying FreeBayes-specific quality filters:"
    log_message "      - QUAL >= 20 (quality score)"
    log_message "      - INFO/DP >= 20 (read depth, no maximum)"
    log_message "      - AF >= 0.05 && AF <= 0.95 (allele frequency)"
    log_message "      - AB >= 0.25 && AB <= 0.75 (allele balance for hets)"
    log_message "      - SAF > 0 && SAR > 0 && SRF > 0 && SRR > 0 (strand support)"
    log_message "      - MQM >= 20 && MQMR >= 20 (mapping quality)"
    
    # Build the filter expression - explicitly use INFO/DP to avoid ambiguity
    local filter_expr="QUAL >= 20 && INFO/DP >= 20"
    
    # Add AF filter if AF field exists
    if bcftools view -h "$input_vcf" | grep -q "##INFO=.*AF"; then
        filter_expr="${filter_expr} && AF >= 0.05 && AF <= 0.95"
    else
        log_message "      - WARNING: AF field not found, skipping AF filter"
    fi
    
    # Add AB filter if AB field exists (allele balance)
    # Check both INFO and FORMAT fields for AB
    local ab_field=""
    if bcftools view -h "$input_vcf" | grep -q "##INFO=.*AB"; then
        ab_field="AB"
    elif bcftools view -h "$input_vcf" | grep -q "##FORMAT=.*AB"; then
        ab_field="FORMAT/AB"
    fi
    
    if [ -n "$ab_field" ]; then
        filter_expr="${filter_expr} && (${ab_field} >= 0.25 && ${ab_field} <= 0.75)"
    else
        log_message "      - WARNING: AB field not found in INFO or FORMAT, skipping allele balance filter"
    fi
    
    # Add strand support filters if available
    local strand_fields_exist=true
    for field in SAF SAR SRF SRR; do
        if ! bcftools view -h "$input_vcf" | grep -q "##INFO=.*${field}"; then
            strand_fields_exist=false
            break
        fi
    done
    
    if [ "$strand_fields_exist" = true ]; then
        filter_expr="${filter_expr} && SAF > 0 && SAR > 0 && SRF > 0 && SRR > 0"
    else
        log_message "      - WARNING: Strand support fields (SAF/SAR/SRF/SRR) not all found, skipping strand filter"
    fi
    
    # Add mapping quality filters if available
    local mq_fields_exist=true
    for field in MQM MQMR; do
        if ! bcftools view -h "$input_vcf" | grep -q "##INFO=.*${field}"; then
            mq_fields_exist=false
            break
        fi
    done
    
    if [ "$mq_fields_exist" = true ]; then
        filter_expr="${filter_expr} && MQM >= 20 && MQMR >= 20"
    else
        log_message "      - WARNING: Mapping quality fields (MQM/MQMR) not found, skipping MQ filter"
    fi
    
    log_message "      - Final filter expression: $filter_expr"
    
    # Apply the filters
    if bcftools filter -i "$filter_expr" -O z -o "$output_vcf" "$input_vcf"; then
        bcftools index -t "$output_vcf"
        
        # Report filtering results
        local input_count=$(bcftools view -H "$input_vcf" | wc -l)
        local output_count=$(bcftools view -H "$output_vcf" | wc -l)
        local filtered_count=$((input_count - output_count))
        
        log_message "      - Input variants: $input_count"
        log_message "      - Passed filters: $output_count"
        log_message "      - Filtered out: $filtered_count"
        
        return 0
    else
        log_message "      - ERROR: Failed to apply FreeBayes filters"
        return 1
    fi
}

# Modified filter_pass_variants function to handle different caller types
filter_pass_variants() {
    local input_vcf="$1"
    local output_vcf="$2"
    
    # Detect caller type
    local caller_type=$(detect_caller_type "$input_vcf")
    log_message "    Detected caller type: $caller_type"
    
    if [[ "$caller_type" == "freebayes" ]]; then
        # Apply FreeBayes-specific quality filters instead of PASS filter
        apply_freebayes_filters "$input_vcf" "$output_vcf"
    else
        # Use PASS filter for other callers (Mutect2, Strelka, etc.)
        log_message "    Applying PASS filter for $caller_type"
        if bcftools view -f PASS "$input_vcf" -O z -o "$output_vcf"; then
            bcftools index -t "$output_vcf"
            
            # Report filtering results
            local input_count=$(bcftools view -H "$input_vcf" | wc -l)
            local output_count=$(bcftools view -H "$output_vcf" | wc -l)
            local filtered_count=$((input_count - output_count))
            
            log_message "      - Input variants: $input_count"
            log_message "      - PASS variants: $output_count"
            log_message "      - Filtered out: $filtered_count"
            
            return 0
        else
            log_message "      - ERROR: Failed to apply PASS filter"
            return 1
        fi
    fi
}

# Function to merge Strelka files if both exist
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

# Modified generate_high_quality_consensus function with enhanced logging
generate_high_quality_consensus() {
    local sample="$1"
    local sample_dir="$2"
    local hq_sample_dir="$3"
    local threshold="$4"
    local sites_file="$5"
    local indexed_vcfs=("${@:6}")
    
    local suffix=""
    if [ "$threshold" -eq 2 ]; then
        suffix="2plus_callers"
    elif [ "$threshold" -eq 3 ]; then
        suffix="3plus_callers"
    else
        suffix="${threshold}plus_callers"
    fi
    
    log_message "  Creating high-quality filtered VCFs for each caller:"
    
    # Create filtered high-quality VCFs
    local pass_vcfs=()
    local pass_count=0
    
    for vcf in "${indexed_vcfs[@]}"; do
        local basename=$(basename "$vcf" .vcf.gz)
        local caller_type=$(detect_caller_type "$vcf")
        local pass_vcf="${sample_dir}/${basename}_HIGH_QUALITY.vcf.gz"
        
        log_message "    Processing $(basename "$vcf") (caller: $caller_type)"
        
        if filter_pass_variants "$vcf" "$pass_vcf"; then
            local variant_count=$(bcftools view -H "$pass_vcf" | wc -l)
            if [ "$variant_count" -gt 0 ]; then
                pass_vcfs+=("$pass_vcf")
                pass_count=$((pass_count + 1))
                log_message "    SUCCESS: Created high-quality VCF: $(basename "$pass_vcf") with $variant_count variants"
            else
                log_message "    WARNING: No variants passed quality filters in $(basename "$vcf")"
                rm -f "$pass_vcf" "${pass_vcf}.tbi"
            fi
        else
            log_message "    ERROR: Failed to create high-quality filtered VCF from $(basename "$vcf")"
        fi
    done
    
    if [ "$pass_count" -lt "$threshold" ]; then
        log_message "  WARNING: Only $pass_count VCFs have high-quality variants, need at least $threshold for consensus"
        # Clean up high-quality VCFs
        for pass_vcf in "${pass_vcfs[@]}"; do
            rm -f "$pass_vcf" "${pass_vcf}.tbi"
        done
        return 1
    fi
    
    log_message "  Running intersection analysis on $pass_count high-quality VCFs..."
    
    # Run bcftools isec on high-quality filtered VCFs
    local hq_isec_dir="${hq_sample_dir}/isec_${suffix}"
    mkdir -p "$hq_isec_dir"
    
    if ! bcftools isec -O v -p "$hq_isec_dir" "${pass_vcfs[@]}"; then
        log_message "  ERROR: bcftools isec failed for high-quality consensus"
        # Clean up
        rm -rf "$hq_isec_dir"
        for pass_vcf in "${pass_vcfs[@]}"; do
            rm -f "$pass_vcf" "${pass_vcf}.tbi"
        done
        return 1
    fi
    
    # Analyze the intersection results for high-quality variants
    local hq_sites_file="${hq_isec_dir}/sites.txt"
    if [ -f "$hq_sites_file" ]; then
        # Count total unique high-quality positions
        local total_hq_positions=$(wc -l < "$hq_sites_file")
        
        # Count variants that appear in at least the threshold number of high-quality files
        local hq_consensus_count=$(awk -v thresh=$threshold '{
            caller_count = 0
            for(i=1; i<=length($5); i++) {
                if(substr($5,i,1) == "1") caller_count++
            }
            if(caller_count >= thresh) print
        }' "$hq_sites_file" | wc -l)
        
        log_message "  High-quality intersection analysis:"
        log_message "    - Total unique high-quality positions: $total_hq_positions"
        log_message "    - Positions in =$threshold callers: $hq_consensus_count"
        
        if [ "$hq_consensus_count" -gt 0 ]; then
            # Extract positions that appear in at least threshold high-quality files
            local hq_consensus_positions="${hq_sample_dir}/hq_consensus_positions_${suffix}.txt"
            awk -v thresh=$threshold '{
                caller_count = 0
                for(i=1; i<=length($5); i++) {
                    if(substr($5,i,1) == "1") caller_count++
                }
                if(caller_count >= thresh) print $1"\t"$2
            }' "$hq_sites_file" > "$hq_consensus_positions"
            
            # Create high-quality consensus VCF using the first high-quality input file as template
            local hq_consensus_vcf="${hq_sample_dir}/high_quality_consensus_${sample}_${suffix}.vcf.gz"
            if bcftools view -R "$hq_consensus_positions" "${pass_vcfs[0]}" -O z -o "$hq_consensus_vcf" 2>/dev/null; then
                bcftools index -t "$hq_consensus_vcf" 2>/dev/null
                
                # Count final variants in high-quality consensus VCF
                local final_hq_variant_count=$(bcftools view -H "$hq_consensus_vcf" | wc -l)
                
                log_message "  SUCCESS: High-quality consensus VCF (=${threshold} high-quality callers) created: $(basename "$hq_consensus_vcf")"
                log_message "  High-quality consensus variants: $final_hq_variant_count"
                
                # Generate stats for high-quality consensus
                log_message "High-quality consensus VCF statistics (=${threshold} high-quality callers):"
                bcftools stats "$hq_consensus_vcf" | grep "^SN" >> "$output_log"
            else
                log_message "  ERROR: Failed to create high-quality consensus VCF for =${threshold} high-quality callers"
                hq_consensus_count=0
            fi
            
            rm -f "$hq_consensus_positions" 2>/dev/null
        else
            log_message "  RESULT: No high-quality consensus variants found for =${threshold} high-quality callers"
        fi
    else
        log_message "  ERROR: No sites.txt file found for high-quality analysis"
        hq_consensus_count=0
    fi
    
    # Clean up temporary files
    rm -rf "$hq_isec_dir"
    for pass_vcf in "${pass_vcfs[@]}"; do
        rm -f "$pass_vcf" "${pass_vcf}.tbi"
    done
    
    return 0
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
    
    # Create sample-specific directories
    sample_dir="${output_base_dir}/${sample}"
    hq_sample_dir="${high_quality_dir}/${sample}"
    mkdir -p "$sample_dir"
    mkdir -p "$hq_sample_dir"
    
    log_message "Processing tumor sample: $sample"
    log_message "Output directory: $sample_dir"
    log_message "High-quality output directory: $hq_sample_dir"
    
    # Merge Strelka files if needed
    merge_strelka_files "$sample" "$sample_dir" vcf_files
    
    log_message "Number of VCF files: ${#vcf_files[@]}"
    
    # Count total variants in input files and show caller-specific filtering info
    total_input_variants=0
    for vcf in "${vcf_files[@]}"; do
        if [ -f "$vcf" ] && [ -s "$vcf" ]; then
            variant_count=$(bcftools view -H "$vcf" | wc -l)
            caller_type=$(detect_caller_type "$vcf")
            
            if [[ "$caller_type" == "freebayes" ]]; then
                log_message "  $(basename "$vcf") ($caller_type): $variant_count total variants (will apply quality filters)"
            else
                pass_count=$(bcftools view -H -f PASS "$vcf" 2>/dev/null | wc -l)
                log_message "  $(basename "$vcf") ($caller_type): $variant_count total variants, $pass_count PASS variants"
            fi
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
        
        # Generate regular consensus files for different thresholds
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
            
            # Generate high-quality consensus (=2 high-quality callers)
            log_message "Generating high-quality consensus (=2 high-quality callers)..."
            generate_high_quality_consensus "$sample" "$sample_dir" "$hq_sample_dir" 2 "$sites_file" "${indexed_vcfs[@]}"
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
            
            # Generate high-quality consensus (=3 high-quality callers)
            log_message "Generating high-quality consensus (=3 high-quality callers)..."
            generate_high_quality_consensus "$sample" "$sample_dir" "$hq_sample_dir" 3 "$sites_file" "${indexed_vcfs[@]}"
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
log_message "High-quality output directory: $high_quality_dir"

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
log_message "Regular consensus output files created:"
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
log_message "High-quality consensus output files created (PASS/quality-filtered variants only):"
for hq_sample_dir in "${high_quality_dir}"/*/; do
    if [ -d "$hq_sample_dir" ]; then
        sample_name=$(basename "$hq_sample_dir")
        log_message "Sample: $sample_name"
        for hq_consensus_file in "$hq_sample_dir"/high_quality_consensus_*.vcf.gz; do
            if [ -f "$hq_consensus_file" ]; then
                file_size=$(ls -lh "$hq_consensus_file" | awk '{print $5}')
                variant_count=$(bcftools view -H "$hq_consensus_file" | wc -l)
                log_message "  $(basename "$hq_consensus_file") (Size: $file_size, Variants: $variant_count)"
            fi
        done
    fi
done

log_message ""
log_message "Processing completed at: $(date)"
log_message "Log file: $output_log"

echo ""
echo "Processing complete! Check the log file: $output_log"
echo "Regular consensus output directory: $output_base_dir"
echo "High-quality consensus output directory: $high_quality_dir"
echo "Total samples processed: $successful_samples"
echo ""
echo "High-quality consensus files contain only variants that:"
echo "  1. Have PASS in the FILTER column (Mutect2/Sarek) OR pass quality filters (FreeBayes)"
echo "  2. Are found in at least 2 (or 3) different variant callers"
echo "  3. Represent the intersection of high-confidence calls"
echo ""
echo "FreeBayes quality filters applied:"
echo "  - QUAL >= 20, DP >= 20, AF 0.05-0.95"
echo "  - AB 0.25-0.75, strand support, MQM/MQMR >= 20"