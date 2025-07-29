#!/bin/bash

module load bcftools

# Set input and output directories
annotated_consensus_dir="/sc/arion/projects/NGSCRC/Work/NRG-GY003/WES/Annotated_Refined_Consensus_VCFs"
annotated_csv="/sc/arion/projects/NGSCRC/Work/Seqera/NRG-GY003_All_Varcalls/csv/annotated_varcalls.csv"
output_dir="/sc/arion/projects/NGSCRC/Work/NRG-GY003/WES/Final_Annotated_Consensus_VCFs"
output_log="${output_dir}/germline_addition_log.txt"

# Germline database cache
germline_db_cache="${output_dir}/germline_database.cache"

# Create output directory
mkdir -p "$output_dir"

echo "Germline Variant Addition Log - $(date)" > "$output_log"
echo "=======================================" >> "$output_log"
echo "" >> "$output_log"

# Function to log messages to both console and file
log_message() {
    echo "$1"
    echo "$1" >> "$output_log"
}

log_message "Starting creation of separate somatic and germline VCF files"
log_message "Annotated consensus directory: $annotated_consensus_dir"
log_message "Annotated VCF CSV: $annotated_csv"
log_message "Output directory: $output_dir"
log_message "Germline database cache: $germline_db_cache"
log_message ""
log_message "Strategy: Create separate VCF files for somatic and germline variants"
log_message "Target callers: FreeBayes, Strelka, HaplotypeCaller"
log_message ""

# Function to check if sample is normal/germline (IMPROVED)
is_normal_sample() {
    local sample_name="$1"
    # Exclude tumor samples (anything with -Tu-)
    if [[ "$sample_name" =~ -Tu- ]]; then
        return 1  # This is a tumor sample
    fi
    return 0  # This is a normal/germline sample
}

# Function to check if caller is one of our target callers
is_target_caller() {
    local caller="$1"
    case "$caller" in
        "freebayes"|"strelka"|"haplotypecaller")
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
    
    if [ ! -f "$cache_file" ]; then
        return 1
    fi
    
    if [ "$csv_file" -nt "$cache_file" ]; then
        return 1
    fi
    
    return 0
}

# Function to save germline database to cache
save_germline_database() {
    local -n db_save_ref=$1
    local cache_file="$2"
    
    log_message "Saving germline database to cache: $cache_file"
    
    > "$cache_file"
    
    for patient in "${!db_save_ref[@]}"; do
        echo "${patient}:${db_save_ref[$patient]}" >> "$cache_file"
    done
    
    log_message "Saved ${#db_save_ref[@]} patient entries to germline cache"
}

# Function to load germline database from cache
load_germline_database() {
    local -n db_load_ref=$1
    local cache_file="$2"
    
    log_message "Loading germline database from cache: $cache_file"
    
    db_load_ref=()
    
    while IFS=':' read -r patient vcf_list; do
        db_load_ref[$patient]="$vcf_list"
    done < "$cache_file"
    
    log_message "Loaded ${#db_load_ref[@]} patient entries from germline cache"
}

# Function to build germline database from CSV
build_germline_database() {
    local -n db_build_ref=$1
    
    log_message "Building germline database from CSV..."
    
    local total_entries=0
    local normal_entries=0
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
        
        # Only process normal samples (exclude tumor samples)
        if ! is_normal_sample "$sample"; then
            continue
        fi
        normal_entries=$((normal_entries + 1))
        
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
            log_message "  Added germline: $patient -> $(basename "$annotatedvcf") [$variantcaller] (sample: $sample)"
        fi
    done < "$annotated_csv"
    
    log_message ""
    log_message "Germline database building summary:"
    log_message "  Total CSV entries: $total_entries"
    log_message "  Normal samples (excluding -Tu-): $normal_entries"
    log_message "  Target callers (freebayes/strelka/haplotypecaller): $target_caller_entries"
    log_message "  Successfully added to database: $added_entries"
    log_message "  Unique patients in database: ${#db_build_ref[@]}"
    log_message ""
}

# Function to get or build germline database
get_germline_database() {
    local -n db_main_ref=$1
    
    if is_cache_valid "$germline_db_cache" "$annotated_csv"; then
        log_message "Using cached germline database..."
        load_germline_database db_main_ref "$germline_db_cache"
        
        # Verify cache integrity
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
                if [ "$total_checked" -ge 5 ]; then
                    break 2
                fi
            done
        done
        
        if [ "$total_checked" -gt 0 ] && [ "$verified_count" -eq "$total_checked" ]; then
            log_message "Germline cache verification passed ($verified_count/$total_checked files exist)"
        else
            log_message "Germline cache verification failed, rebuilding..."
            build_germline_database db_main_ref
            save_germline_database db_main_ref "$germline_db_cache"
        fi
    else
        log_message "Germline cache is stale or missing, building new database..."
        build_germline_database db_main_ref
        save_germline_database db_main_ref "$germline_db_cache"
    fi
}

# Function to extract pathogenic germline variants (OPTIMIZED - targeted frequency filtering with 5% threshold)
extract_pathogenic_germline() {
    local patient_id="$1"
    local output_germline_vcf="$2"
    shift 2
    local vcf_files=("$@")
    
    log_message "      Extracting pathogenic germline variants from ${#vcf_files[@]} normal VCFs"
    log_message "      Requiring variants to be present in >=2 callers for consensus"
    log_message "      Enhanced filtering: HIGH impact variants in protein-coding genes only"
    log_message "      Population frequency filter: >5% (from VEP CSQ fields 42,48,58,69)"
    log_message "      Normalizing variants before consensus calling"
    
    local temp_dir="/tmp/germline_${patient_id}_$$"
    mkdir -p "$temp_dir"
    
    # Reference genome path
    local reference_genome="/sc/arion/projects/NGSCRC/Resources/gatk_hg38/gatk_bundle_2024/Homo_sapiens_assembly38.fasta"
    
    # FIXED: Better module loading for bgzip
    local can_normalize=true
    
    # Try to load modules that provide bgzip
    if ! command -v bgzip >/dev/null 2>&1; then
        log_message "        bgzip not found, trying to load modules..."
        
        # Try different module names that might provide bgzip
        for module_name in "htslib" "samtools" "bcftools" "HTSlib" "SAMtools"; do
            if module load "$module_name" 2>/dev/null; then
                log_message "        Loaded module: $module_name"
                if command -v bgzip >/dev/null 2>&1; then
                    log_message "        ? bgzip now available"
                    break
                fi
            fi
        done
    fi
    
    # Final check for bgzip and reference genome
    if ! command -v bgzip >/dev/null 2>&1; then
        log_message "        WARNING: bgzip still not available after trying modules"
        log_message "        Available modules: $(module avail 2>&1 | grep -i hts || echo 'none found')"
        can_normalize=false
    elif [ ! -f "$reference_genome" ]; then
        log_message "        WARNING: Reference genome not found at $reference_genome"
        can_normalize=false
    else
        log_message "        ? bgzip and reference genome available"
    fi
    
    local filtered_vcfs=()
    local normalized_vcfs=()
    local total_pathogenic=0
    local caller_names=()
    
    # Process each germline VCF
    for vcf in "${vcf_files[@]}"; do
        local caller_name="unknown"
        local annotated_csv="/sc/arion/projects/NGSCRC/Work/Seqera/NRG-GY003_All_Varcalls/csv/annotated_varcalls.csv"
        while IFS=',' read -r patient sample variantcaller annotatedvcf; do
            if [[ "$patient" == "patient" ]]; then
                continue
            fi
    
            patient=$(echo "$patient" | tr -d '"' | tr -d '\r' | xargs)
            variantcaller=$(echo "$variantcaller" | tr -d '"' | tr -d '\r' | xargs)
            annotatedvcf=$(echo "$annotatedvcf" | tr -d '"' | tr -d '\r' | xargs)
    
            if [[ "$annotatedvcf" == "$vcf" ]]; then
                caller_name="$variantcaller"
                break
            fi
        done < "$annotated_csv"
        local filtered_vcf="${temp_dir}/pathogenic_${caller_name}.vcf"
        local normalized_vcf="${temp_dir}/normalized_${caller_name}.vcf"
        
        log_message "        Processing $(basename "$vcf") [$caller_name]"
        
        # Ensure VCF is indexed
        if [ ! -f "${vcf}.tbi" ]; then
            bcftools index -t "$vcf"
        fi
        
        # Test what annotation fields are actually available
        local has_ann=0
        local has_csq=0
        local has_clnsig=0
        
        if bcftools view -h "$vcf" | grep -q "##INFO=.*ANN"; then
            has_ann=1
        fi
        
        if bcftools view -h "$vcf" | grep -q "##INFO=.*CSQ"; then
            has_csq=1
        fi
        
        if bcftools view -h "$vcf" | grep -q "##INFO=.*CLNSIG"; then
            has_clnsig=1
        fi
        
        # Build filter based on what's available
        local filter_parts=()
        
        # Add HIGH impact filters if annotation fields exist
        if [ "$has_ann" -eq 1 ]; then
            filter_parts+=('(INFO/ANN ~ "HIGH")')
            filter_parts+=('(INFO/ANN ~ "stop_gained")')
            filter_parts+=('(INFO/ANN ~ "frameshift")')
            filter_parts+=('(INFO/ANN ~ "splice_acceptor")')
            filter_parts+=('(INFO/ANN ~ "splice_donor")')
            filter_parts+=('(INFO/ANN ~ "start_lost")')
            filter_parts+=('(INFO/ANN ~ "stop_lost")')
            # Add cancer genes with precise matching
            filter_parts+=('(INFO/ANN ~ "\\|BRCA1\\|")')
            filter_parts+=('(INFO/ANN ~ "\\|BRCA2\\|")')
            filter_parts+=('(INFO/ANN ~ "\\|TP53\\|")')
            filter_parts+=('(INFO/ANN ~ "\\|PTEN\\|")')
            filter_parts+=('(INFO/ANN ~ "\\|ATM\\|")')
            filter_parts+=('(INFO/ANN ~ "\\|CHEK2\\|")')
            filter_parts+=('(INFO/ANN ~ "\\|PALB2\\|")')
            filter_parts+=('(INFO/ANN ~ "\\|APC\\|")')
            filter_parts+=('(INFO/ANN ~ "\\|MUTYH\\|")')
        fi
        
        if [ "$has_csq" -eq 1 ]; then
            filter_parts+=('(INFO/CSQ ~ "HIGH")')
            filter_parts+=('(INFO/CSQ ~ "stop_gained")')
            filter_parts+=('(INFO/CSQ ~ "frameshift_variant")')
            filter_parts+=('(INFO/CSQ ~ "splice_acceptor_variant")')
            filter_parts+=('(INFO/CSQ ~ "splice_donor_variant")')
            filter_parts+=('(INFO/CSQ ~ "start_lost")')
            filter_parts+=('(INFO/CSQ ~ "stop_lost")')
            # Add cancer genes with precise matching
            filter_parts+=('(INFO/CSQ ~ "\\|BRCA1\\|")')
            filter_parts+=('(INFO/CSQ ~ "\\|BRCA2\\|")')
            filter_parts+=('(INFO/CSQ ~ "\\|TP53\\|")')
            filter_parts+=('(INFO/CSQ ~ "\\|PTEN\\|")')
            filter_parts+=('(INFO/CSQ ~ "\\|ATM\\|")')
            filter_parts+=('(INFO/CSQ ~ "\\|CHEK2\\|")')
            filter_parts+=('(INFO/CSQ ~ "\\|PALB2\\|")')
            filter_parts+=('(INFO/CSQ ~ "\\|APC\\|")')
            filter_parts+=('(INFO/CSQ ~ "\\|MUTYH\\|")')
        fi
        
        if [ "$has_clnsig" -eq 1 ]; then
            filter_parts+=('(INFO/CLNSIG ~ "Pathogenic")')
            filter_parts+=('(INFO/CLNSIG ~ "Likely_pathogenic")')
            filter_parts+=('(INFO/CLNSIG ~ "pathogenic")')
            filter_parts+=('(INFO/CLNSIG ~ "likely_pathogenic")')
        fi
        
        # If no annotation fields available, skip this VCF
        if [ ${#filter_parts[@]} -eq 0 ]; then
            log_message "          No annotation fields available, skipping"
            continue
        fi
        
        # Join filter parts with OR
        local pathogenic_filter=""
        for i in "${!filter_parts[@]}"; do
            if [ $i -eq 0 ]; then
                pathogenic_filter="${filter_parts[i]}"
            else
                pathogenic_filter="${pathogenic_filter} || ${filter_parts[i]}"
            fi
        done
        
        # Wrap in parentheses
        pathogenic_filter="($pathogenic_filter)"
        
        # Apply broad filter first
        if bcftools filter -i "$pathogenic_filter" -o "${temp_dir}/temp_${caller_name}.vcf" "$vcf" 2>/dev/null; then
            
            # Apply optimized AWK filter with targeted frequency filtering (5% threshold)
            awk -v has_clnsig="$has_clnsig" '
            BEGIN { FS="\t"; OFS="\t" }
            /^#/ { print; next }
            {
                keep = 0
                
                # Always keep cancer gene variants
                if ($8 ~ /\|BRCA1\|| \|BRCA2\|| \|TP53\|| \|PTEN\|| \|ATM\|| \|CHEK2\|| \|PALB2\|| \|APC\|| \|MUTYH\|/) {
                    keep = 1
                }
                
                # Keep ClinVar pathogenic variants only if CLNSIG field exists
                if (has_clnsig == 1 && $8 ~ /CLNSIG.*[Pp]athogenic/ && $8 !~ /benign/) {
                    keep = 1
                }
                
                # For HIGH impact variants, check benign status and population frequency
                if ($8 ~ /HIGH/) {
                    is_benign = 0
                    if ($8 ~ /benign/) is_benign = 1
                    
                    # OPTIMIZED: Check specific frequency fields from CSQ (5% threshold)
                    is_common = 0
                    
                    if ($8 ~ /CSQ=/) {
                        csq_field = $8
                        gsub(/.*CSQ=/, "", csq_field)
                        gsub(/;.*/, "", csq_field)
                        
                        # Split CSQ annotations by comma
                        split(csq_field, csq_annotations, ",")
                        for (csq_i in csq_annotations) {
                            split(csq_annotations[csq_i], csq_fields, "|")
                            
                            # Check key population frequency fields (targeted approach)
                            # Field 42: AF (1KG), Field 48: gnomADe_AF, Field 58: gnomADg_AF, Field 69: MAX_AF
                            key_freq_fields[1] = csq_fields[42]  # AF (1000 Genomes)
                            key_freq_fields[2] = csq_fields[48]  # gnomADe_AF  
                            key_freq_fields[3] = csq_fields[58]  # gnomADg_AF
                            key_freq_fields[4] = csq_fields[69]  # MAX_AF
                            
                            for (freq_idx in key_freq_fields) {
                                freq_val = key_freq_fields[freq_idx]
                                if (freq_val != "" && freq_val ~ /^[0-9]*\.?[0-9]+$/ && freq_val > 0.05) {  # 5% threshold
                                    is_common = 1
                                    break
                                }
                            }
                            if (is_common) break
                        }
                    }
                    
                    if (!is_benign && !is_common) {
                        if ($8 ~ /ANN=/) {
                            ann_field = $8
                            gsub(/.*ANN=/, "", ann_field)
                            gsub(/;.*/, "", ann_field)
                            
                            split(ann_field, annotations, ",")
                            for (j in annotations) {
                                split(annotations[j], fields, "|")
                                if (fields[3] == "HIGH" && fields[8] == "protein_coding") {
                                    keep = 1
                                    break
                                }
                            }
                        }
                        
                        if ($8 ~ /CSQ=/) {
                            csq_field = $8
                            gsub(/.*CSQ=/, "", csq_field)
                            gsub(/;.*/, "", csq_field)
                            
                            split(csq_field, annotations, ",")
                            for (j in annotations) {
                                split(annotations[j], fields, "|")
                                if (fields[3] == "HIGH" && fields[8] == "protein_coding") {
                                    keep = 1
                                    break
                                }
                            }
                        }
                    }
                }
                
                # For specific variant types (similar optimized frequency check)
                if ($8 ~ /stop_gained|frameshift|splice_acceptor|splice_donor|start_lost|stop_lost/) {
                    is_benign = 0
                    if ($8 ~ /benign/) is_benign = 1
                    
                    # OPTIMIZED: Check specific frequency fields from CSQ (5% threshold)
                    is_common = 0
                    
                    if ($8 ~ /CSQ=/) {
                        csq_field = $8
                        gsub(/.*CSQ=/, "", csq_field)
                        gsub(/;.*/, "", csq_field)
                        
                        split(csq_field, csq_annotations, ",")
                        for (csq_i in csq_annotations) {
                            split(csq_annotations[csq_i], csq_fields, "|")
                            
                            # Check key population frequency fields
                            key_freq_fields[1] = csq_fields[42]  # AF (1000 Genomes)
                            key_freq_fields[2] = csq_fields[48]  # gnomADe_AF  
                            key_freq_fields[3] = csq_fields[58]  # gnomADg_AF
                            key_freq_fields[4] = csq_fields[69]  # MAX_AF
                            
                            for (freq_idx in key_freq_fields) {
                                freq_val = key_freq_fields[freq_idx]
                                if (freq_val != "" && freq_val ~ /^[0-9]*\.?[0-9]+$/ && freq_val > 0.05) {  # 5% threshold
                                    is_common = 1
                                    break
                                }
                            }
                            if (is_common) break
                        }
                    }
                    
                    if (!is_benign && !is_common) {
                        if ($8 ~ /ANN=/) {
                            ann_field = $8
                            gsub(/.*ANN=/, "", ann_field)
                            gsub(/;.*/, "", ann_field)
                            
                            split(ann_field, annotations, ",")
                            for (j in annotations) {
                                split(annotations[j], fields, "|")
                                if ((fields[2] ~ /stop_gained|frameshift|splice_acceptor|splice_donor|start_lost|stop_lost/) && fields[8] == "protein_coding") {
                                    keep = 1
                                    break
                                }
                            }
                        }
                    }
                }
                
                if (keep) print
            }' "${temp_dir}/temp_${caller_name}.vcf" > "$filtered_vcf"
            
            local variant_count=$(grep -v "^#" "$filtered_vcf" | wc -l)
            if [ "$variant_count" -gt 0 ]; then
                log_message "          Found $variant_count pathogenic variants"
                
                # Normalize variants if possible
                if [ "$can_normalize" = true ]; then
                    log_message "          Normalizing variants..."
                    
                    # ALTERNATIVE: Use gzip instead of bgzip if bgzip still fails
                    local compressed_vcf="${filtered_vcf}.gz"
                    if bgzip -c "$filtered_vcf" > "$compressed_vcf" 2>/dev/null; then
                        # bgzip worked
                        :
                    elif gzip -c "$filtered_vcf" > "$compressed_vcf" 2>/dev/null; then
                        # fallback to gzip
                        log_message "          Using gzip instead of bgzip"
                    else
                        log_message "          Compression failed, skipping normalization"
                        cp "$filtered_vcf" "$normalized_vcf"
                        normalized_vcfs+=("$normalized_vcf")
                        caller_names+=("$caller_name")
                        total_pathogenic=$((total_pathogenic + variant_count))
                        continue
                    fi
                    
                    if bcftools index "$compressed_vcf" 2>/dev/null; then
                        # Normalize: left-align and split multiallelic sites
                        if bcftools norm -f "$reference_genome" -m -any "$compressed_vcf" > "$normalized_vcf" 2>/dev/null; then
                            local normalized_count=$(grep -v "^#" "$normalized_vcf" | wc -l)
                            log_message "          Normalized to $normalized_count variants"
                            
                            normalized_vcfs+=("$normalized_vcf")
                            caller_names+=("$caller_name")
                            total_pathogenic=$((total_pathogenic + normalized_count))
                        else
                            log_message "          Normalization failed, using original"
                            cp "$filtered_vcf" "$normalized_vcf"
                            normalized_vcfs+=("$normalized_vcf")
                            caller_names+=("$caller_name")
                            total_pathogenic=$((total_pathogenic + variant_count))
                        fi
                        
                        # Clean up temp files
                        rm -f "$compressed_vcf" "${compressed_vcf}.csi" "${compressed_vcf}.tbi"
                    else
                        log_message "          Indexing failed, using original"
                        cp "$filtered_vcf" "$normalized_vcf"
                        normalized_vcfs+=("$normalized_vcf")
                        caller_names+=("$caller_name")
                        total_pathogenic=$((total_pathogenic + variant_count))
                        rm -f "$compressed_vcf"
                    fi
                else
                    # If no normalization possible, use original
                    cp "$filtered_vcf" "$normalized_vcf"
                    normalized_vcfs+=("$normalized_vcf")
                    caller_names+=("$caller_name")
                    total_pathogenic=$((total_pathogenic + variant_count))
                fi
            else
                log_message "          No pathogenic variants found"
            fi
            
            # Clean up temp file
            rm -f "${temp_dir}/temp_${caller_name}.vcf"
        else
            log_message "          ERROR: Filter command failed"
        fi
    done
    
    # Apply multi-caller consensus filtering using NORMALIZED variants
    if [ ${#normalized_vcfs[@]} -lt 2 ]; then
        log_message "      Only ${#normalized_vcfs[@]} caller(s) with pathogenic variants - insufficient for consensus"
        log_message "      Creating empty germline VCF (requires >=2 callers)"
        
        # Create empty VCF
        if [ ${#vcf_files[@]} -gt 0 ]; then
            bcftools view -h "${vcf_files[0]}" | bcftools view -O z -o "$output_germline_vcf"
            bcftools index -t "$output_germline_vcf"
        fi
        
        # Clean up
        rm -rf "$temp_dir"
        return 0
    fi
    
    log_message "        Applying multi-caller consensus on normalized variants (>=2 callers required)"
    
    # Use bcftools isec to find variants present in >=2 callers (using normalized variants)
    local isec_dir="${temp_dir}/isec_consensus"
    mkdir -p "$isec_dir"
    
    # Convert normalized VCFs to compressed format for isec
    local compressed_vcfs=()
    for i in "${!normalized_vcfs[@]}"; do
        local norm_vcf="${normalized_vcfs[i]}"
        local caller="${caller_names[i]}"
        local compressed_vcf="${isec_dir}/norm_${caller}.vcf.gz"
        
        # Try bgzip first, then gzip
        if bgzip -c "$norm_vcf" > "$compressed_vcf" 2>/dev/null || gzip -c "$norm_vcf" > "$compressed_vcf" 2>/dev/null; then
            if bcftools index "$compressed_vcf" 2>/dev/null; then
                compressed_vcfs+=("$compressed_vcf")
            fi
        fi
    done
    
    # Run isec to find overlapping variants
    if [ ${#compressed_vcfs[@]} -ge 2 ] && bcftools isec -p "$isec_dir" -n +2 "${compressed_vcfs[@]}" 2>/dev/null; then
        
        # Collect consensus variants from all intersection files
        local consensus_vcfs=()
        local consensus_variants=0
        
        for i in $(seq 0 $((${#compressed_vcfs[@]} - 1))); do
            local isec_file="${isec_dir}/000${i}.vcf"
            if [ -f "$isec_file" ] && [ -s "$isec_file" ]; then
                # Convert to compressed and index
                bcftools view -O z -o "${isec_file}.gz" "$isec_file"
                bcftools index -t "${isec_file}.gz"
                consensus_vcfs+=("${isec_file}.gz")
                
                local variant_count=$(bcftools view -H "${isec_file}.gz" | wc -l)
                consensus_variants=$((consensus_variants + variant_count))
                log_message "          ${caller_names[i]}: $variant_count consensus variants"
            fi
        done
        
        if [ ${#consensus_vcfs[@]} -gt 0 ]; then
            # Use the first consensus source
            cp "${consensus_vcfs[0]}" "$output_germline_vcf"
            cp "${consensus_vcfs[0]}.tbi" "${output_germline_vcf}.tbi"
            
            # Count final consensus variants
            local final_consensus_count=$(bcftools view -H "$output_germline_vcf" | wc -l)
            
            log_message "      Enhanced multi-caller consensus results:"
            log_message "        - Total pathogenic variants found: $total_pathogenic"
            log_message "        - Consensus variants (>=2 callers): $final_consensus_count"
            log_message "        - Callers contributing: ${#normalized_vcfs[@]}"
            if [ "$can_normalize" = true ]; then
                log_message "        - Applied: Normalization + protein-coding filter + population frequency filter (5%)"
            else
                log_message "        - Applied: protein-coding filter + population frequency filter (5%) (no normalization)"
            fi
            
        else
            log_message "        No consensus variants found (no overlap between callers)"
            # Create empty VCF
            bcftools view -h "${vcf_files[0]}" | bcftools view -O z -o "$output_germline_vcf"
            bcftools index -t "$output_germline_vcf"
        fi
        
    else
        log_message "        ERROR: bcftools isec failed for consensus analysis"
        # Fallback: use first normalized VCF if available
        if [ ${#normalized_vcfs[@]} -gt 0 ]; then
            bcftools view -O z -o "$output_germline_vcf" "${normalized_vcfs[0]}"
            bcftools index -t "$output_germline_vcf"
        fi
    fi
    
    # Clean up
    rm -rf "$temp_dir"
    return 0
}

# Function to process a single sample (UPDATED - one germline VCF per patient)
process_sample_with_separate_files() {
    local sample_dir="$1"
    local patient_id=$(basename "$sample_dir")
    local -n germline_db_ref=$2
    
    log_message "Processing sample with separate file creation: $patient_id"
    
    # Create output directory
    local output_sample_dir="${output_dir}/${patient_id}"
    mkdir -p "$output_sample_dir"
    
    # Check if we have germline VCFs for this patient
    if [[ -z "${germline_db_ref[$patient_id]}" ]]; then
        log_message "  WARNING: No germline VCFs found for patient $patient_id"
        log_message "  Creating somatic-only VCFs"
        
        # Create somatic-only VCFs
        for consensus_vcf in "$sample_dir"/annotated_*.vcf.gz; do
            if [ -f "$consensus_vcf" ]; then
                local basename=$(basename "$consensus_vcf" .vcf.gz)
                local somatic_vcf="${output_sample_dir}/final_${basename}_SOMATIC.vcf.gz"
                
                # Tag as somatic-only
                bcftools view -h "$consensus_vcf" > "${output_sample_dir}/temp_somatic.vcf"
                bcftools view -H "$consensus_vcf" | awk 'BEGIN { FS=OFS="\t" } {
                    if ($3 == ".") {
                        $3 = "SOMATIC"
                    } else {
                        $3 = $3 "_SOMATIC"
                    }
                    print
                }' >> "${output_sample_dir}/temp_somatic.vcf"
                
                bcftools view -O z -o "$somatic_vcf" "${output_sample_dir}/temp_somatic.vcf"
                bcftools index -t "$somatic_vcf"
                rm "${output_sample_dir}/temp_somatic.vcf"
                
                local variant_count=$(bcftools view -H "$somatic_vcf" | wc -l)
                log_message "    Created: $(basename "$somatic_vcf") ($variant_count somatic variants)"
            fi
        done
        return 1
    fi
    
    local germline_vcf_list="${germline_db_ref[$patient_id]}"
    local germline_vcfs=($germline_vcf_list)
    
    log_message "  Found ${#germline_vcfs[@]} germline VCFs for patient $patient_id"
    
    # Extract pathogenic germline variants ONCE for this patient
    local pathogenic_germline="${output_sample_dir}/temp_consensus_pathogenic_germline_${patient_id}.vcf.gz"
    extract_pathogenic_germline "$patient_id" "$pathogenic_germline" "${germline_vcfs[@]}"
    
    local processed_files=0
    local successful_files=0
    local germline_variant_count=$(bcftools view -H "$pathogenic_germline" | wc -l)
    
    # Create somatic VCFs for each consensus threshold
    for consensus_vcf in "$sample_dir"/annotated_*.vcf.gz; do
        if [ -f "$consensus_vcf" ]; then
            processed_files=$((processed_files + 1))
            
            local basename=$(basename "$consensus_vcf" .vcf.gz)
            local somatic_vcf="${output_sample_dir}/final_${basename}_SOMATIC.vcf.gz"
            
            log_message "    Creating somatic VCF: $(basename "$somatic_vcf")"
            
            # Create somatic VCF
            bcftools view -h "$consensus_vcf" > "${output_sample_dir}/temp_somatic.vcf"
            bcftools view -H "$consensus_vcf" | awk 'BEGIN { FS=OFS="\t" } {
                if ($3 == ".") {
                    $3 = "SOMATIC"
                } else {
                    $3 = $3 "_SOMATIC"
                }
                print
            }' >> "${output_sample_dir}/temp_somatic.vcf"
            
            bcftools view -O z -o "$somatic_vcf" "${output_sample_dir}/temp_somatic.vcf"
            bcftools index -t "$somatic_vcf"
            rm "${output_sample_dir}/temp_somatic.vcf"
            
            local somatic_count=$(bcftools view -H "$somatic_vcf" | wc -l)
            successful_files=$((successful_files + 1))
            
            log_message "      Created: $(basename "$somatic_vcf") ($somatic_count somatic variants)"
        fi
    done
    
    # Create one germline VCF (if variants exist)
    if [ "$germline_variant_count" -gt 0 ]; then
        local germline_output="${output_sample_dir}/final_consensus_pathogenic_germline_${patient_id}.vcf.gz"
        
        log_message "    Creating consensus germline VCF: $(basename "$germline_output")"
        
        # Tag germline variants
        bcftools view -h "$pathogenic_germline" > "${output_sample_dir}/temp_germline.vcf"
        bcftools view -H "$pathogenic_germline" | awk 'BEGIN { FS=OFS="\t" } {
            if ($3 == ".") {
                $3 = "GERMLINE_PATHOGENIC"
            } else {
                $3 = $3 "_GERMLINE_PATHOGENIC"
            }
            print
        }' >> "${output_sample_dir}/temp_germline.vcf"
        
        bcftools view -O z -o "$germline_output" "${output_sample_dir}/temp_germline.vcf"
        bcftools index -t "$germline_output"
        rm "${output_sample_dir}/temp_germline.vcf"
        
        log_message "      Created: $(basename "$germline_output") ($germline_variant_count consensus germline variants)"
    else
        log_message "    No consensus germline variants found (requires >=2 callers)"
    fi
    
    # Clean up intermediate files
    rm -f "$pathogenic_germline"
    rm -f "${pathogenic_germline}.tbi"
    
    log_message "  Sample summary: $successful_files/$processed_files consensus VCFs processed successfully"
    log_message "---"
    
    return 0
}

# Main processing
log_message "Scanning for annotated consensus VCF directories..."

if [ ! -d "$annotated_consensus_dir" ]; then
    log_message "ERROR: Annotated consensus directory $annotated_consensus_dir does not exist!"
    exit 1
fi

if [ ! -f "$annotated_csv" ]; then
    log_message "ERROR: Annotated VCF CSV file $annotated_csv does not exist!"
    exit 1
fi

# Get germline database
declare -A germline_database
get_germline_database germline_database

# Process samples
total_samples=0
successful_samples=0

for sample_dir in "$annotated_consensus_dir"/*/; do
    if [ -d "$sample_dir" ]; then
        total_samples=$((total_samples + 1))
        
        if process_sample_with_separate_files "$sample_dir" germline_database; then
            successful_samples=$((successful_samples + 1))
        fi
    fi
done

# Count total files created (CORRECTED)
somatic_vcfs=0
germline_vcfs=0
total_somatic_variants=0
total_germline_variants=0

for output_sample_dir in "$output_dir"/*/; do
    if [ -d "$output_sample_dir" ]; then
        for somatic_vcf in "$output_sample_dir"/final_*_SOMATIC.vcf.gz; do
            if [ -f "$somatic_vcf" ]; then
                somatic_vcfs=$((somatic_vcfs + 1))
                somatic_count=$(bcftools view -H "$somatic_vcf" | wc -l)
                total_somatic_variants=$((total_somatic_variants + somatic_count))
            fi
        done
        
        # UPDATED: Look for the new germline file pattern
        for germline_vcf in "$output_sample_dir"/final_consensus_pathogenic_germline_*.vcf.gz; do
            if [ -f "$germline_vcf" ]; then
                germline_vcfs=$((germline_vcfs + 1))
                germline_count=$(bcftools view -H "$germline_vcf" | wc -l)
                total_germline_variants=$((total_germline_variants + germline_count))
            fi
        done
    fi
done

# Final summary
log_message ""
log_message "======================================="
log_message "SEPARATE FILE CREATION COMPLETE"
log_message "======================================="
log_message "Total samples found: $total_samples"
log_message "Successfully processed samples: $successful_samples"
log_message "Failed samples: $((total_samples - successful_samples))"
log_message "Total somatic VCF files created: $somatic_vcfs"
log_message "Total germline VCF files created: $germline_vcfs"
log_message "Total somatic variants: $total_somatic_variants"
log_message "Total pathogenic germline variants: $total_germline_variants"
log_message ""
log_message "File naming convention:"
log_message "  - Somatic VCFs: final_*_SOMATIC.vcf.gz"
log_message "  - Germline VCFs: final_consensus_pathogenic_germline_PATIENT.vcf.gz"
log_message "  - Variants tagged in ID field for easy identification"
log_message ""
log_message "Target callers used: FreeBayes, Strelka, HaplotypeCaller"
log_message "Pathogenic germline criteria used:"
log_message "  - HIGH impact variants (SnpEff/VEP)"
log_message "  - ClinVar Pathogenic/Likely pathogenic"
log_message "  - Protein-truncating variants (stop-gained, frameshift, splice)"
log_message "  - Cancer predisposition genes (BRCA1/2, TP53, PTEN, ATM, CHEK2, PALB2, APC, MUTYH)"
log_message "  - Multi-caller consensus (>=2 callers required)"
log_message "  - Population frequency filter: >5% (from VEP CSQ fields 42,48,58,69)"
log_message "  - Excludes benign variants"
log_message ""
log_message "Created VCF files:"

for output_sample_dir in "$output_dir"/*/; do
    if [ -d "$output_sample_dir" ]; then
        sample_name=$(basename "$output_sample_dir")
        log_message "Sample: $sample_name"
        
        for vcf_file in "$output_sample_dir"/final_*.vcf.gz; do
            if [ -f "$vcf_file" ]; then
                file_size=$(ls -lh "$vcf_file" | awk '{print $5}')
                variant_count=$(bcftools view -H "$vcf_file" | wc -l)
                
                if [[ "$vcf_file" =~ _SOMATIC\.vcf\.gz$ ]]; then
                    log_message "  $(basename "$vcf_file") (Size: $file_size, Somatic variants: $variant_count)"
                elif [[ "$vcf_file" =~ final_consensus_pathogenic_germline_.*\.vcf\.gz$ ]]; then
                    log_message "  $(basename "$vcf_file") (Size: $file_size, Germline variants: $variant_count)"
                fi
            fi
        done
    fi
done

log_message ""
log_message "Processing completed at: $(date)"
log_message "Log file: $output_log"

echo ""
echo "Separate file creation complete!"
echo "Input directory: $annotated_consensus_dir"
echo "Output directory: $output_dir"
echo "Log file: $output_log"
echo ""
echo "? Created separate somatic and germline VCF files"
echo "? Somatic VCFs: $somatic_vcfs files with $total_somatic_variants variants"
echo "? Germline VCFs: $germline_vcfs files with $total_germline_variants variants"
echo "? Ready for separate MAF conversion"
echo ""
echo "File naming:"
echo "  - Somatic: final_*_SOMATIC.vcf.gz"
echo "  - Germline: final_consensus_pathogenic_germline_PATIENT.vcf.gz"
echo ""
echo "For MAF conversion, process each file type separately:"
echo "  1. Convert somatic VCFs ? somatic MAF"
echo "  2. Convert germline VCFs ? germline MAF"
echo "  3. Combine MAFs if needed"