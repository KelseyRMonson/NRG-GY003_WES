#!/bin/bash

module load bcftools

# Set input and output directories
input_dir="/sc/arion/projects/NGSCRC/Work/NRG-GY003/WES/High_Quality_Consensus_VCFs"
output_dir="/sc/arion/projects/NGSCRC/Work/NRG-GY003/WES/Refined_Consensus_VCFs"
output_log="${output_dir}/refined_vcf_processing_log.txt"

# Create output directory
mkdir -p "$output_dir"

echo "Refined VCF Processing Log - $(date)" > "$output_log"
echo "=======================================" >> "$output_log"
echo "" >> "$output_log"

# Function to log messages to both console and file
log_message() {
    echo "$1"
    echo "$1" >> "$output_log"
}

# Define quality filters - SIMPLIFIED
MIN_TUMOR_DEPTH=20     # Minimum tumor sample depth (FORMAT/DP) - increased from 10 to 20
MIN_AF=0.05           # Minimum tumor allele frequency (strict)
MAX_AF=0.95           # Maximum tumor allele frequency

log_message "Starting refined quality filtering of high-quality consensus VCFs"
log_message "Input directory: $input_dir"
log_message "Output directory: $output_dir"
log_message ""
log_message "Quality filters (VCFs already filtered to PASS variants):"
log_message "  - Tumor FORMAT/DP >= $MIN_TUMOR_DEPTH (tumor sample depth)"
log_message "  - Tumor FORMAT/AF: $MIN_AF <= AF <= $MAX_AF (tumor allele frequency - strict)"
log_message ""

# Function to analyze VCF structure (simplified)
analyze_vcf_structure() {
    local vcf_file="$1"
    
    log_message "    Analyzing VCF structure for $(basename "$vcf_file"):"
    
    # Show sample names
    local samples=$(bcftools query -l "$vcf_file" | tr '\n' ' ')
    log_message "      Samples: $samples"
    
    # Count total variants
    local total_variants=$(bcftools view -H "$vcf_file" | wc -l)
    log_message "      Total variants: $total_variants"
    
    # Identify tumor sample
    local sample_list=($(bcftools query -l "$vcf_file"))
    local tumor_idx=1  # Default to second sample (index 1)
    
    for i in "${!sample_list[@]}"; do
        if [[ "${sample_list[i]}" =~ Tu ]]; then
            tumor_idx=$i
            break
        fi
    done
    
    log_message "      Using sample ${sample_list[tumor_idx]} (index $tumor_idx) as tumor sample"
}

# Function to apply refined quality filters with simplified breakdown
apply_refined_filters() {
    local input_vcf="$1"
    local output_vcf="$2"
    local sample_name="$3"
    
    log_message "  Processing: $(basename "$input_vcf")"
    
    # Analyze VCF structure first
    analyze_vcf_structure "$input_vcf"
    
    # Get sample names to identify tumor sample
    local sample_list=($(bcftools query -l "$input_vcf"))
    local tumor_idx=1  # Default to second sample
    
    for i in "${!sample_list[@]}"; do
        if [[ "${sample_list[i]}" =~ Tu ]]; then
            tumor_idx=$i
            break
        fi
    done
    
    log_message "    Applying tumor-focused filters:"
    
    # Create AWK script for filtering with simplified breakdown
    local filter_script="/tmp/vcf_filter_${sample_name}_$$.awk"
    
    cat > "$filter_script" << EOF
BEGIN { 
    FS="\t"; OFS="\t"
    MIN_TUMOR_DEPTH=$MIN_TUMOR_DEPTH
    MIN_AF=$MIN_AF
    MAX_AF=$MAX_AF
    TUMOR_IDX=$tumor_idx
    
    # Counters for breakdown
    total_variants=0
    passed_all=0
    failed_tumor_dp=0
    failed_tumor_af_low=0
    failed_tumor_af_high=0
}

# Print header lines
/^#/ { print; next }

# Process variant lines
{
    total_variants++
    
    # Parse FORMAT field to find AF and DP positions
    split(\$9, format_fields, ":")
    af_pos = 0
    dp_pos = 0
    for(i=1; i<=length(format_fields); i++) {
        if(format_fields[i] == "AF") af_pos = i
        if(format_fields[i] == "DP") dp_pos = i
    }
    
    # Get tumor sample column (10 + tumor_idx)
    tumor_col = 10 + TUMOR_IDX
    if(tumor_col > NF) tumor_col = NF  # Fallback to last column
    
    # Parse tumor sample data
    split(\$tumor_col, tumor_data, ":")
    
    tumor_af = (af_pos > 0 && af_pos <= length(tumor_data)) ? tumor_data[af_pos] : "."
    tumor_dp = (dp_pos > 0 && dp_pos <= length(tumor_data)) ? tumor_data[dp_pos] : "."
    
    # Check each filter criterion
    pass_tumor_dp = (tumor_dp != "." && tumor_dp >= MIN_TUMOR_DEPTH)
    pass_tumor_af = (tumor_af != "." && tumor_af >= MIN_AF && tumor_af <= MAX_AF)
    
    # Count failures for each criterion
    if(!pass_tumor_dp) failed_tumor_dp++
    if(tumor_af != "." && tumor_af < MIN_AF) failed_tumor_af_low++
    if(tumor_af != "." && tumor_af > MAX_AF) failed_tumor_af_high++
    
    # Apply filters - only tumor-based criteria
    if(pass_tumor_dp && pass_tumor_af) {
        print
        passed_all++
    }
}

END {
    printf "# FILTERING BREAKDOWN:\n" > "/dev/stderr"
    printf "# Total variants: %d\n", total_variants > "/dev/stderr"
    printf "# Passed all filters: %d (%.1f%%)\n", passed_all, (passed_all/total_variants)*100 > "/dev/stderr"
    printf "# Failed filters: %d (%.1f%%)\n", total_variants-passed_all, ((total_variants-passed_all)/total_variants)*100 > "/dev/stderr"
    printf "#\n" > "/dev/stderr"
    printf "# Failure breakdown by criterion:\n" > "/dev/stderr"
    printf "#   Tumor DP < %d: %d variants\n", MIN_TUMOR_DEPTH, failed_tumor_dp > "/dev/stderr"
    printf "#   Tumor AF < %.2f: %d variants\n", MIN_AF, failed_tumor_af_low > "/dev/stderr"
    printf "#   Tumor AF > %.2f: %d variants\n", MAX_AF, failed_tumor_af_high > "/dev/stderr"
}
EOF
    
    # Count input variants
    local input_count=$(bcftools view -H "$input_vcf" | wc -l)
    
    # Apply the AWK filter and capture the breakdown
    local filter_output=$(bcftools view "$input_vcf" | awk -f "$filter_script" 2>&1 >/tmp/filtered_output_$$.vcf)
    
    # Convert to compressed VCF
    if bcftools view -O z -o "$output_vcf" /tmp/filtered_output_$$.vcf; then
        bcftools index -t "$output_vcf"
        
        # Count output variants
        local output_count=$(bcftools view -H "$output_vcf" | wc -l)
        local filtered_count=$((input_count - output_count))
        
        if [ "$input_count" -gt 0 ]; then
            local retention_rate=$(awk "BEGIN {printf \"%.2f\", $output_count/$input_count*100}")
        else
            local retention_rate="0.00"
        fi
        
        log_message "    Results:"
        log_message "      - Input variants: $input_count"
        log_message "      - Passed tumor-focused filters: $output_count"
        log_message "      - Filtered out: $filtered_count"
        log_message "      - Retention rate: ${retention_rate}%"
        log_message ""
        log_message "    Filtering breakdown by criterion:"
        
        # Parse and log the breakdown
        echo "$filter_output" | grep "^#" | sed 's/^# /      /' | while read line; do
            log_message "$line"
        done
        
        # Clean up
        rm -f "$filter_script" /tmp/filtered_output_$$.vcf
        
        return 0
    else
        log_message "    ERROR: Failed to apply refined filters"
        rm -f "$filter_script" /tmp/filtered_output_$$.vcf
        return 1
    fi
}

# Function to process a single sample directory
process_sample() {
    local sample_dir="$1"
    local sample_name=$(basename "$sample_dir")
    
    log_message "Processing sample: $sample_name"
    
    # Create output directory for this sample
    local output_sample_dir="${output_dir}/${sample_name}"
    mkdir -p "$output_sample_dir"
    
    local processed_files=0
    local successful_files=0
    
    # Process all high-quality consensus VCFs in the sample directory
    for vcf_file in "$sample_dir"/high_quality_consensus_*.vcf.gz; do
        if [ -f "$vcf_file" ]; then
            processed_files=$((processed_files + 1))
            
            # Generate output filename
            local basename=$(basename "$vcf_file" .vcf.gz)
            local output_vcf="${output_sample_dir}/refined_${basename}.vcf.gz"
            
            if apply_refined_filters "$vcf_file" "$output_vcf" "$sample_name"; then
                successful_files=$((successful_files + 1))
            fi
        fi
    done
    
    if [ "$processed_files" -eq 0 ]; then
        log_message "  WARNING: No high-quality consensus VCF files found in $sample_dir"
        return 1
    fi
    
    log_message "  Sample summary: $successful_files/$processed_files files processed successfully"
    log_message "---"
    
    return 0
}

# Main processing
log_message "Scanning for sample directories..."

if [ ! -d "$input_dir" ]; then
    log_message "ERROR: Input directory $input_dir does not exist!"
    exit 1
fi

# Initialize counters
total_samples=0
successful_samples=0
successful_vcfs=0

# Process each sample directory
for sample_dir in "$input_dir"/*/; do
    if [ -d "$sample_dir" ]; then
        total_samples=$((total_samples + 1))
        
        if process_sample "$sample_dir"; then
            successful_samples=$((successful_samples + 1))
        fi
    fi
done

# Count total VCFs processed
for output_sample_dir in "$output_dir"/*/; do
    if [ -d "$output_sample_dir" ]; then
        for refined_vcf in "$output_sample_dir"/refined_*.vcf.gz; do
            if [ -f "$refined_vcf" ]; then
                successful_vcfs=$((successful_vcfs + 1))
            fi
        done
    fi
done

# Final summary
log_message ""
log_message "======================================="
log_message "FINAL SUMMARY"
log_message "======================================="
log_message "Total samples found: $total_samples"
log_message "Successfully processed samples: $successful_samples"
log_message "Failed samples: $((total_samples - successful_samples))"
log_message "Total refined VCF files created: $successful_vcfs"
log_message ""
log_message "Refined VCF files created with tumor-focused quality filters:"

for output_sample_dir in "$output_dir"/*/; do
    if [ -d "$output_sample_dir" ]; then
        sample_name=$(basename "$output_sample_dir")
        log_message "Sample: $sample_name"
        
        for refined_vcf in "$output_sample_dir"/refined_*.vcf.gz; do
            if [ -f "$refined_vcf" ]; then
                file_size=$(ls -lh "$refined_vcf" | awk '{print $5}')
                variant_count=$(bcftools view -H "$refined_vcf" | wc -l)
                log_message "  $(basename "$refined_vcf") (Size: $file_size, Variants: $variant_count)"
            fi
        done
    fi
done

log_message ""
log_message "Tumor-focused quality filters applied:"
log_message "  - Tumor FORMAT/DP >= $MIN_TUMOR_DEPTH (tumor sample depth)"
log_message "  - Tumor FORMAT/AF: $MIN_AF <= AF <= $MAX_AF (tumor allele frequency - strict)"
log_message ""
log_message "Note: Removed total depth filter - focusing only on tumor sample metrics"
log_message ""
log_message "Processing completed at: $(date)"
log_message "Log file: $output_log"

echo ""
echo "Refined quality filtering complete!"
echo "Input directory: $input_dir"
echo "Output directory: $output_dir"
echo "Log file: $output_log"
echo "Total refined VCF files created: $successful_vcfs"
echo ""
echo "Applied tumor-focused filters: DP >= 20 and AF 0.05-0.95 in tumor sample only."