#!/bin/bash

module load bcftools

# Set input and output directories
refined_vcf_dir="/sc/arion/projects/NGSCRC/Work/NRG-GY003/WES/Refined_Consensus_VCFs"
annotated_csv="/sc/arion/projects/NGSCRC/Work/Seqera/NRG-GY003_All_Varcalls/csv/annotated_varcalls.csv"
output_dir="/sc/arion/projects/NGSCRC/Work/NRG-GY003/WES/Annotated_Refined_Consensus_VCFs"
output_log="${output_dir}/annotation_processing_log.txt"

# Database cache file
annotation_db_cache="${output_dir}/annotation_database.cache"

# Create output directory
mkdir -p "$output_dir"

echo "Annotation Processing Log - $(date)" > "$output_log"
echo "=======================================" >> "$output_log"
echo "" >> "$output_log"

# Function to log messages to both console and file
log_message() {
    echo "$1"
    echo "$1" >> "$output_log"
}

log_message "Starting FAST annotation transfer using bcftools isec"
log_message "Refined VCF directory: $refined_vcf_dir"
log_message "Annotated VCF CSV: $annotated_csv"
log_message "Output directory: $output_dir"
log_message "Database cache: $annotation_db_cache"
log_message ""

# Function to check if sample is tumor-related
is_tumor_sample() {
    local sample_name="$1"
    if [[ "$sample_name" =~ -Tu- ]]; then
        return 0
    fi
    return 1
}

# Function to check if caller is one of our target callers
is_target_caller() {
    local caller="$1"
    case "$caller" in
        "freebayes"|"mutect2"|"strelka")
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

# Function to check if cache is valid
is_cache_valid() {
    local cache_file="$1"
    local csv_file="$2"
    
    # Check if cache exists
    if [ ! -f "$cache_file" ]; then
        return 1  # Cache doesn't exist
    fi
    
    # Check if CSV is newer than cache
    if [ "$csv_file" -nt "$cache_file" ]; then
        return 1  # CSV is newer, cache is stale
    fi
    
    return 0  # Cache is valid
}

# Function to save annotation database to cache
save_annotation_database() {
    local -n db_save_ref=$1
    local cache_file="$2"
    
    log_message "Saving annotation database to cache: $cache_file"
    
    # Create cache file with patient -> VCF mappings
    > "$cache_file"  # Clear file
    
    for patient in "${!db_save_ref[@]}"; do
        echo "${patient}:${db_save_ref[$patient]}" >> "$cache_file"
    done
    
    log_message "Saved ${#db_save_ref[@]} patient entries to cache"
}

# Function to load annotation database from cache
load_annotation_database() {
    local -n db_load_ref=$1
    local cache_file="$2"
    
    log_message "Loading annotation database from cache: $cache_file"
    
    # Clear existing database
    db_load_ref=()
    
    # Load from cache
    while IFS=':' read -r patient vcf_list; do
        db_load_ref[$patient]="$vcf_list"
    done < "$cache_file"
    
    log_message "Loaded ${#db_load_ref[@]} patient entries from cache"
}

# Function to build annotation database from CSV (streamlined)
build_annotation_database() {
    local -n db_build_ref=$1
    
    log_message "Building streamlined annotation database from CSV..."
    
    local total_entries=0
    local tumor_entries=0
    local target_caller_entries=0
    local added_entries=0
    
    while IFS=',' read -r patient sample variantcaller annotatedvcf; do
        if [[ "$patient" == "patient" ]]; then
            continue
        fi
        
        total_entries=$((total_entries + 1))
        
        patient=$(echo "$patient" | tr -d '"' | tr -d '\r' | xargs)
        sample=$(echo "$sample" | tr -d '"' | tr -d '\r' | xargs)
        variantcaller=$(echo "$variantcaller" | tr -d '"' | tr -d '\r' | xargs)
        annotatedvcf=$(echo "$annotatedvcf" | tr -d '"' | tr -d '\r' | xargs)
        
        if ! is_tumor_sample "$sample"; then
            continue
        fi
        tumor_entries=$((tumor_entries + 1))
        
        if ! is_target_caller "$variantcaller"; then
            continue
        fi
        target_caller_entries=$((target_caller_entries + 1))
        
        if [[ "$annotatedvcf" =~ snpEff_VEP\.ann\.vcf\.gz$ ]] && [[ -f "$annotatedvcf" ]]; then
            if [[ -n "${db_build_ref[$patient]}" ]]; then
                db_build_ref[$patient]="${db_build_ref[$patient]} $annotatedvcf"
            else
                db_build_ref[$patient]="$annotatedvcf"
            fi
            added_entries=$((added_entries + 1))
            log_message "  Added: $patient -> $(basename "$annotatedvcf") [$variantcaller]"
        fi
    done < "$annotated_csv"
    
    log_message ""
    log_message "Database building summary:"
    log_message "  Total CSV entries: $total_entries"
    log_message "  Tumor samples: $tumor_entries"
    log_message "  Target callers (freebayes/mutect2/strelka): $target_caller_entries"
    log_message "  Successfully added to database: $added_entries"
    log_message "  Unique patients in database: ${#db_build_ref[@]}"
    log_message ""
}

# Function to get or build annotation database
get_annotation_database() {
    local -n db_main_ref=$1
    
    if is_cache_valid "$annotation_db_cache" "$annotated_csv"; then
        log_message "Using cached annotation database..."
        load_annotation_database db_main_ref "$annotation_db_cache"
        
        # Verify cache integrity by checking a few files
        local verified_count=0
        local total_checked=0
        for patient in "${!db_main_ref[@]}"; do
            local vcf_list="${db_main_ref[$patient]}"
            local vcfs=($vcf_list)
            for vcf in "${vcfs[@]}"; do
                total_checked=$((total_checked + 1))
                if [ -f "$vcf" ]; then
                    verified_count=$((verified_count + 1))
                fi
                # Only check first few files for speed
                if [ "$total_checked" -ge 10 ]; then
                    break 2
                fi
            done
        done
        
        if [ "$total_checked" -gt 0 ] && [ "$verified_count" -eq "$total_checked" ]; then
            log_message "Cache verification passed ($verified_count/$total_checked files exist)"
        else
            log_message "Cache verification failed ($verified_count/$total_checked files exist), rebuilding..."
            build_annotation_database db_main_ref
            save_annotation_database db_main_ref "$annotation_db_cache"
        fi
    else
        log_message "Cache is stale or missing, building new annotation database..."
        build_annotation_database db_main_ref
        save_annotation_database db_main_ref "$annotation_db_cache"
    fi
}

# Function to fast annotate using bcftools isec (FIXED)
fast_annotate_consensus_vcf() {
    local consensus_vcf="$1"
    local patient_id="$2"
    local annotated_vcf_list="$3"
    local output_vcf="$4"
    
    log_message "    FAST annotating $(basename "$consensus_vcf") for patient $patient_id"
    
    local annotated_vcfs=($annotated_vcf_list)
    log_message "      Using ${#annotated_vcfs[@]} tumor-specific annotated VCFs"
    
    local temp_dir="/tmp/fast_annotation_${patient_id}_$$"
    mkdir -p "$temp_dir"
    
    # Step 1: Use bcftools isec to find intersections quickly
    local isec_dir="${temp_dir}/isec"
    mkdir -p "$isec_dir"
    
    # Prepare all VCFs for intersection (consensus + annotated)
    local all_vcfs=("$consensus_vcf")
    for vcf in "${annotated_vcfs[@]}"; do
        if [ ! -f "${vcf}.tbi" ]; then
            bcftools index -t "$vcf"
        fi
        all_vcfs+=("$vcf")
    done
    
    # Ensure consensus VCF is indexed
    if [ ! -f "${consensus_vcf}.tbi" ]; then
        bcftools index -t "$consensus_vcf"
    fi
    
    log_message "      Running bcftools isec for fast intersection..."
    
    # Run isec to find overlapping positions
    if bcftools isec -p "$isec_dir" -n +2 "${all_vcfs[@]}" 2>/dev/null; then
        
        # Find which files have overlaps with consensus (file 0)
        local annotation_sources=()
        local total_matches=0
        
        for i in $(seq 1 $((${#annotated_vcfs[@]}))); do
            local isec_file="${isec_dir}/000${i}.vcf"
            if [ -f "$isec_file" ] && [ -s "$isec_file" ]; then
                # Convert to compressed and index
                bcftools view -O z -o "${isec_file}.gz" "$isec_file"
                bcftools index -t "${isec_file}.gz"
                annotation_sources+=("${isec_file}.gz")
                local match_count=$(bcftools view -H "${isec_file}.gz" | wc -l)
                total_matches=$((total_matches + match_count))
                log_message "        Found $match_count matches in $(basename "${annotated_vcfs[$((i-1))]}")"
            fi
        done
        
        if [ ${#annotation_sources[@]} -gt 0 ]; then
            # Handle merging with duplicate sample names
            local merged_annotations="${temp_dir}/merged_annotations.vcf.gz"
            
            if [ ${#annotation_sources[@]} -eq 1 ]; then
                cp "${annotation_sources[0]}" "$merged_annotations"
                cp "${annotation_sources[0]}.tbi" "${merged_annotations}.tbi"
                log_message "      Using single annotation source"
            else
                log_message "      Merging ${#annotation_sources[@]} annotation sources (handling duplicate sample names)"
                
                # Use --force-samples to handle duplicate sample names and ensure proper compression
                if bcftools merge --force-samples -m none -O z -o "$merged_annotations" "${annotation_sources[@]}" 2>/dev/null; then
                    # Verify the file is properly compressed
                    if bcftools view -H "$merged_annotations" >/dev/null 2>&1; then
                        bcftools index -t "$merged_annotations"
                        log_message "        Successfully merged annotation sources"
                    else
                        log_message "        Merge produced invalid file, using first source only"
                        cp "${annotation_sources[0]}" "$merged_annotations"
                        cp "${annotation_sources[0]}.tbi" "${merged_annotations}.tbi"
                    fi
                else
                    log_message "        Merge failed, using first annotation source only"
                    cp "${annotation_sources[0]}" "$merged_annotations"
                    cp "${annotation_sources[0]}.tbi" "${merged_annotations}.tbi"
                fi
            fi
            
            # Extract annotation headers
            local annotation_headers="${temp_dir}/annotation_headers.txt"
            if bcftools view -h "$merged_annotations" | grep "^##INFO" | grep -E "(ANN|CSQ|EFF|snpEff|VEP)" > "$annotation_headers"; then
                local header_count=$(wc -l < "$annotation_headers")
                log_message "        Found $header_count annotation field definitions"
            else
                log_message "        No annotation headers found"
                echo "" > "$annotation_headers"
            fi
            
            # Get annotation fields
            local annotation_fields=""
            if [ -s "$annotation_headers" ]; then
                annotation_fields=$(grep "^##INFO" "$annotation_headers" | \
                                   sed 's/^##INFO=<ID=\([^,]*\).*/\1/' | \
                                   tr '\n' ',' | sed 's/,$//')
            fi
            
            if [ -n "$annotation_fields" ]; then
                log_message "      Transferring fields: $(echo "$annotation_fields" | tr ',' ' ')"
                
                # Annotate consensus VCF
                if bcftools annotate \
                    -a "$merged_annotations" \
                    -c "$annotation_fields" \
                    -h "$annotation_headers" \
                    -O z -o "$output_vcf" \
                    "$consensus_vcf" 2>/dev/null; then
                    
                    bcftools index -t "$output_vcf"
                    
                    # Count results
                    local input_count=$(bcftools view -H "$consensus_vcf" | wc -l)
                    local output_count=$(bcftools view -H "$output_vcf" | wc -l)
                    local annotated_count=$(bcftools view -H "$output_vcf" | grep -E "(ANN=|CSQ=)" | wc -l)
                    local annotation_rate=$(awk "BEGIN {printf \"%.1f\", $annotated_count/$input_count*100}")
                    
                    log_message "      FAST Results:"
                    log_message "        - Input variants: $input_count"
                    log_message "        - Output variants: $output_count"
                    log_message "        - Annotated variants: $annotated_count ($annotation_rate%)"
                    log_message "        - Total matches found: $total_matches"
                    
                else
                    log_message "      ERROR: Failed to annotate"
                    cp "$consensus_vcf" "$output_vcf"
                    bcftools index -t "$output_vcf"
                fi
            else
                log_message "      WARNING: No annotation fields found"
                cp "$consensus_vcf" "$output_vcf"
                bcftools index -t "$output_vcf"
            fi
        else
            log_message "      No overlapping variants found"
            cp "$consensus_vcf" "$output_vcf"
            bcftools index -t "$output_vcf"
        fi
    else
        log_message "      ERROR: bcftools isec failed"
        cp "$consensus_vcf" "$output_vcf"
        bcftools index -t "$output_vcf"
    fi
    
    # Clean up
    rm -rf "$temp_dir"
    return 0
}

# Function to process a single sample directory
process_sample() {
    local sample_dir="$1"
    local patient_id=$(basename "$sample_dir")
    local -n db_process_ref=$2
    
    log_message "Processing sample: $patient_id"
    
    if [[ -z "${db_process_ref[$patient_id]}" ]]; then
        log_message "  WARNING: No tumor-specific annotated VCFs found for patient $patient_id"
        return 1
    fi
    
    local annotated_vcf_list="${db_process_ref[$patient_id]}"
    local annotated_vcfs=($annotated_vcf_list)
    
    log_message "  Found ${#annotated_vcfs[@]} tumor-specific annotated VCFs for patient $patient_id"
    
    local output_sample_dir="${output_dir}/${patient_id}"
    mkdir -p "$output_sample_dir"
    
    local processed_files=0
    local successful_files=0
    
    for refined_vcf in "$sample_dir"/refined_*.vcf.gz; do
        if [ -f "$refined_vcf" ]; then
            processed_files=$((processed_files + 1))
            
            local basename=$(basename "$refined_vcf" .vcf.gz)
            local output_vcf="${output_sample_dir}/annotated_${basename}.vcf.gz"
            
            if fast_annotate_consensus_vcf "$refined_vcf" "$patient_id" "$annotated_vcf_list" "$output_vcf"; then
                successful_files=$((successful_files + 1))
            fi
        fi
    done
    
    if [ "$processed_files" -eq 0 ]; then
        log_message "  WARNING: No refined consensus VCF files found in $sample_dir"
        return 1
    fi
    
    log_message "  Sample summary: $successful_files/$processed_files files processed successfully"
    log_message "---"
    
    return 0
}

# Main processing
log_message "Scanning for refined consensus VCF directories..."

if [ ! -d "$refined_vcf_dir" ]; then
    log_message "ERROR: Refined VCF directory $refined_vcf_dir does not exist!"
    exit 1
fi

if [ ! -f "$annotated_csv" ]; then
    log_message "ERROR: Annotated VCF CSV file $annotated_csv does not exist!"
    exit 1
fi

# Get annotation database (from cache or build new)
declare -A annotation_database
get_annotation_database annotation_database

# Process samples
total_samples=0
successful_samples=0

for sample_dir in "$refined_vcf_dir"/*/; do
    if [ -d "$sample_dir" ]; then
        total_samples=$((total_samples + 1))
        
        if process_sample "$sample_dir" annotation_database; then
            successful_samples=$((successful_samples + 1))
        fi
    fi
done

# Count total annotated VCFs created and gather statistics
successful_vcfs=0
total_input_variants=0
total_annotated_variants=0

for output_sample_dir in "$output_dir"/*/; do
    if [ -d "$output_sample_dir" ]; then
        for annotated_vcf in "$output_sample_dir"/annotated_*.vcf.gz; do
            if [ -f "$annotated_vcf" ]; then
                successful_vcfs=$((successful_vcfs + 1))
                variant_count=$(bcftools view -H "$annotated_vcf" | wc -l)
                annotation_count=$(bcftools view -H "$annotated_vcf" | grep -E "(ANN=|CSQ=)" | wc -l)
                total_input_variants=$((total_input_variants + variant_count))
                total_annotated_variants=$((total_annotated_variants + annotation_count))
            fi
        done
    fi
done

# Calculate overall annotation rate
overall_annotation_rate="0.0"
if [ "$total_input_variants" -gt 0 ]; then
    overall_annotation_rate=$(awk "BEGIN {printf \"%.1f\", $total_annotated_variants/$total_input_variants*100}")
fi

# Final summary
log_message ""
log_message "======================================="
log_message "FAST PROCESSING COMPLETE"
log_message "======================================="
log_message "Total samples found: $total_samples"
log_message "Successfully processed samples: $successful_samples"
log_message "Failed samples: $((total_samples - successful_samples))"
log_message "Total annotated VCF files created: $successful_vcfs"
log_message "Total variants processed: $total_input_variants"
log_message "Total variants with annotations: $total_annotated_variants"
log_message "Overall annotation rate: ${overall_annotation_rate}%"
log_message ""
log_message "FAST processing used bcftools isec for bulk intersection"
log_message "Focused on tumor samples with FreeBayes, Mutect2, and Strelka only"
log_message "Database cache saved to: $annotation_db_cache"
log_message ""

echo ""
echo "FAST annotation transfer complete!"
echo "Database cache saved to: $annotation_db_cache"
echo "Next run will be faster using cached database!"
echo "Total annotated consensus VCF files created: $successful_vcfs"
echo "Overall annotation rate: ${overall_annotation_rate}%"