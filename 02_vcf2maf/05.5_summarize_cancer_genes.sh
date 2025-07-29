#!/bin/bash

# Script to summarize cancer gene variants per patient from germline VCFs (FINAL - handles missing frequencies)
module load bcftools

output_dir="/sc/arion/projects/NGSCRC/Work/NRG-GY003/WES/Final_Annotated_Consensus_VCFs"
summary_file="${output_dir}/cancer_gene_germline_summary.txt"
detailed_file="${output_dir}/cancer_gene_germline_detailed.csv"

# Cancer genes to analyze
cancer_genes=("BRCA1" "BRCA2" "TP53" "PTEN" "ATM" "CHEK2" "PALB2" "APC" "MUTYH")

echo "=== Cancer Gene Germline Variant Summary ===" > "$summary_file"
echo "Date: $(date)" >> "$summary_file"
echo "Source: Consensus pathogenic germline VCFs (>=2 callers)" >> "$summary_file"
echo "Filters applied: HIGH impact + protein-coding + <5% population frequency + not benign" >> "$summary_file"
echo "" >> "$summary_file"

# Create detailed CSV header
echo "Patient,Gene,Chromosome,Position,Ref,Alt,Consequence,Impact,HGVSc,HGVSp,Protein_Change,Frequency_1KG,Frequency_gnomADe,Frequency_gnomADg,Max_Frequency,Clinical_Significance,Frequency_Status,Callers_Supporting" > "$detailed_file"

echo "Analyzing cancer gene variants in germline VCFs..."
echo ""

# Initialize counters
declare -A gene_patient_counts
declare -A gene_variant_counts
declare -A consequence_counts
declare -A frequency_status_counts
total_patients_with_variants=0
total_patients_analyzed=0

# Process each patient directory
for patient_dir in "$output_dir"/*/; do
    if [ -d "$patient_dir" ]; then
        patient_id=$(basename "$patient_dir")
        total_patients_analyzed=$((total_patients_analyzed + 1))
        
        # Look for germline VCF
        germline_vcf="${patient_dir}/final_consensus_pathogenic_germline_${patient_id}.vcf.gz"
        
        if [ -f "$germline_vcf" ]; then
            echo "Processing patient: $patient_id"
            
            # Check if patient has any cancer gene variants
            patient_has_cancer_variants=false
            patient_summary=""
            
            # Analyze each cancer gene
            for gene in "${cancer_genes[@]}"; do
                # Count variants for this gene using more precise matching
                variant_count=0
                if bcftools view -H "$germline_vcf" | grep -q "|${gene}|"; then
                    variant_count=$(bcftools view -H "$germline_vcf" | grep -c "|${gene}|")
                fi
                
                if [ "$variant_count" -gt 0 ]; then
                    patient_has_cancer_variants=true
                    gene_variant_counts[$gene]=$((${gene_variant_counts[$gene]:-0} + variant_count))
                    
                    # Add patient to gene count if not already counted
                    if [[ ! "${gene_patient_counts[$gene]}" =~ $patient_id ]]; then
                        gene_patient_counts[$gene]="${gene_patient_counts[$gene]:-} $patient_id"
                    fi
                    
                    patient_summary="${patient_summary}  $gene: $variant_count variants\n"
                    
                    # Extract detailed information for each variant
                    bcftools view -H "$germline_vcf" | grep "|${gene}|" | while read line; do
                        chrom=$(echo "$line" | cut -f1)
                        pos=$(echo "$line" | cut -f2)
                        ref=$(echo "$line" | cut -f4)
                        alt=$(echo "$line" | cut -f5)
                        info=$(echo "$line" | cut -f8)
                        
                        # Initialize variables
                        consequence="unknown"
                        impact="unknown"
                        hgvsc="unknown"
                        hgvsp="unknown"
                        protein_change="unknown"
                        freq_1kg="not_available"
                        freq_gnomade="not_available"
                        freq_gnomadg="not_available"
                        max_freq="not_available"
                        clinical_sig="unknown"
                        frequency_status="no_population_data"
                        
                        # Extract from ANN field first (SnpEff)
                        if [[ "$info" =~ ANN=([^;]+) ]]; then
                            ann="${BASH_REMATCH[1]}"
                            # Split by comma and find the annotation for this gene
                            IFS=',' read -ra ann_parts <<< "$ann"
                            for ann_part in "${ann_parts[@]}"; do
                                if [[ "$ann_part" =~ \|${gene}\| ]]; then
                                    IFS='|' read -ra fields <<< "$ann_part"
                                    consequence="${fields[1]:-unknown}"
                                    impact="${fields[2]:-unknown}"
                                    hgvsc="${fields[9]:-unknown}"
                                    hgvsp="${fields[10]:-unknown}"
                                    protein_change="${fields[15]:-unknown}"
                                    break
                                fi
                            done
                        fi
                        
                        # Extract frequency information from CSQ field (VEP)
                        # Try both gene-specific and first annotation
                        if [[ "$info" =~ CSQ=([^;]+) ]]; then
                            csq="${BASH_REMATCH[1]}"
                            IFS=',' read -ra csq_parts <<< "$csq"
                            
                            # First try gene-specific annotation
                            for csq_part in "${csq_parts[@]}"; do
                                if [[ "$csq_part" =~ \|${gene}\| ]]; then
                                    IFS='|' read -ra csq_fields <<< "$csq_part"
                                    
                                    # Check if this annotation has frequency data
                                    if [ -n "${csq_fields[41]}" ] || [ -n "${csq_fields[47]}" ] || [ -n "${csq_fields[57]}" ] || [ -n "${csq_fields[68]}" ]; then
                                        freq_1kg="${csq_fields[41]:-not_available}"
                                        freq_gnomade="${csq_fields[47]:-not_available}"
                                        freq_gnomadg="${csq_fields[57]:-not_available}"
                                        max_freq="${csq_fields[68]:-not_available}"
                                        frequency_status="gene_specific_data"
                                    fi
                                    
                                    clinical_sig="${csq_fields[71]:-unknown}"
                                    
                                    # Also get consequence/impact from CSQ if ANN didn't work
                                    if [ "$consequence" = "unknown" ]; then
                                        consequence="${csq_fields[1]:-unknown}"
                                        impact="${csq_fields[2]:-unknown}"
                                        hgvsc="${csq_fields[10]:-unknown}"
                                        hgvsp="${csq_fields[11]:-unknown}"
                                    fi
                                    break
                                fi
                            done
                            
                            # If no gene-specific frequency data, try first annotation (canonical)
                            if [ "$frequency_status" = "no_population_data" ] && [ ${#csq_parts[@]} -gt 0 ]; then
                                IFS='|' read -ra csq_fields <<< "${csq_parts[0]}"
                                if [ -n "${csq_fields[41]}" ] || [ -n "${csq_fields[47]}" ] || [ -n "${csq_fields[57]}" ] || [ -n "${csq_fields[68]}" ]; then
                                    freq_1kg="${csq_fields[41]:-not_available}"
                                    freq_gnomade="${csq_fields[47]:-not_available}"
                                    freq_gnomadg="${csq_fields[57]:-not_available}"
                                    max_freq="${csq_fields[68]:-not_available}"
                                    frequency_status="canonical_annotation_data"
                                fi
                            fi
                        fi
                        
                        # Clean up frequency values
                        [ -z "$freq_1kg" ] && freq_1kg="not_available"
                        [ -z "$freq_gnomade" ] && freq_gnomade="not_available"
                        [ -z "$freq_gnomadg" ] && freq_gnomadg="not_available"
                        [ -z "$max_freq" ] && max_freq="not_available"
                        [ -z "$clinical_sig" ] && clinical_sig="unknown"
                        
                        # Determine frequency status for summary
                        if [[ "$freq_1kg" != "not_available" ]] || [[ "$freq_gnomade" != "not_available" ]] || [[ "$freq_gnomadg" != "not_available" ]] || [[ "$max_freq" != "not_available" ]]; then
                            frequency_status="has_population_data"
                            # Check if it's rare (all frequencies < 0.01 or not_available)
                            is_rare=true
                            for freq_val in "$freq_1kg" "$freq_gnomade" "$freq_gnomadg" "$max_freq"; do
                                if [[ "$freq_val" != "not_available" ]] && [[ $(echo "$freq_val > 0.01" | bc -l 2>/dev/null) == "1" ]]; then
                                    is_rare=false
                                    break
                                fi
                            done
                            if [ "$is_rare" = true ]; then
                                frequency_status="rare_variant"
                            else
                                frequency_status="common_variant"
                            fi
                        else
                            frequency_status="no_population_data"
                        fi
                        
                        # Count frequency status
                        frequency_status_counts["$frequency_status"]=$((${frequency_status_counts["$frequency_status"]:-0} + 1))
                        
                        # Count consequences for summary
                        if [[ "$consequence" =~ stop_gained|frameshift|splice_acceptor|splice_donor|start_lost|stop_lost ]]; then
                            consequence_counts["$consequence"]=$((${consequence_counts["$consequence"]:-0} + 1))
                        fi
                        
                        # Determine supporting callers (simplified)
                        callers_supporting=">=2_callers"
                        
                        # Write to detailed CSV (escape commas in fields)
                        consequence_clean=$(echo "$consequence" | tr ',' ';')
                        hgvsc_clean=$(echo "$hgvsc" | tr ',' ';')
                        hgvsp_clean=$(echo "$hgvsp" | tr ',' ';')
                        protein_change_clean=$(echo "$protein_change" | tr ',' ';')
                        clinical_sig_clean=$(echo "$clinical_sig" | tr ',' ';')
                        
                        echo "$patient_id,$gene,$chrom,$pos,$ref,$alt,$consequence_clean,$impact,$hgvsc_clean,$hgvsp_clean,$protein_change_clean,$freq_1kg,$freq_gnomade,$freq_gnomadg,$max_freq,$clinical_sig_clean,$frequency_status,$callers_supporting" >> "$detailed_file"
                    done
                fi
            done
            
            if [ "$patient_has_cancer_variants" = true ]; then
                total_patients_with_variants=$((total_patients_with_variants + 1))
                echo "Patient $patient_id:" >> "$summary_file"
                echo -e "$patient_summary" >> "$summary_file"
            fi
        else
            echo "  No germline VCF found for patient $patient_id"
        fi
    fi
done

# Create summary statistics
echo "" >> "$summary_file"
echo "=== SUMMARY STATISTICS ===" >> "$summary_file"
echo "Total patients analyzed: $total_patients_analyzed" >> "$summary_file"
echo "Patients with cancer gene variants: $total_patients_with_variants" >> "$summary_file"
if [ "$total_patients_analyzed" -gt 0 ]; then
    percentage=$(( (total_patients_with_variants * 100) / total_patients_analyzed ))
    echo "Percentage with cancer gene variants: $percentage%" >> "$summary_file"
fi
echo "" >> "$summary_file"

echo "Summary by gene:" >> "$summary_file"
echo "================" >> "$summary_file"
for gene in "${cancer_genes[@]}"; do
    variant_count=${gene_variant_counts[$gene]:-0}
    patient_list="${gene_patient_counts[$gene]:-}"
    patient_count=$(echo "$patient_list" | wc -w)
    
    if [ "$patient_count" -gt 0 ]; then
        percentage=$(( (patient_count * 100) / total_patients_analyzed ))
        echo "$gene: $variant_count variants in $patient_count patients ($percentage%)" >> "$summary_file"
        echo "  Patients with $gene variants:$patient_list" >> "$summary_file"
    else
        echo "$gene: No variants found" >> "$summary_file"
    fi
    echo "" >> "$summary_file"
done

# Add frequency status summary
echo "=== FREQUENCY DATA AVAILABILITY ===" >> "$summary_file"
echo "Population frequency data status:" >> "$summary_file"
for status in "no_population_data" "rare_variant" "common_variant"; do
    count=${frequency_status_counts[$status]:-0}
    if [ "$count" -gt 0 ]; then
        echo "  $status: $count variants" >> "$summary_file"
    fi
done
echo "" >> "$summary_file"
echo "Note: 'no_population_data' indicates very rare/novel variants not in major databases" >> "$summary_file"
echo "This is common and expected for pathogenic cancer predisposition variants" >> "$summary_file"
echo "" >> "$summary_file"

# Create a patient-by-gene matrix
echo "=== PATIENT-GENE MATRIX ===" >> "$summary_file"
printf "%-12s" "Patient" >> "$summary_file"
for gene in "${cancer_genes[@]}"; do
    printf "%-8s" "$gene" >> "$summary_file"
done
printf "%-8s\n" "Total" >> "$summary_file"

echo "================================================================" >> "$summary_file"

for patient_dir in "$output_dir"/*/; do
    if [ -d "$patient_dir" ]; then
        patient_id=$(basename "$patient_dir")
        germline_vcf="${patient_dir}/final_consensus_pathogenic_germline_${patient_id}.vcf.gz"
        
        if [ -f "$germline_vcf" ]; then
            patient_total=0
            gene_counts=()
            
            for gene in "${cancer_genes[@]}"; do
                variant_count=0
                if bcftools view -H "$germline_vcf" | grep -q "|${gene}|"; then
                    variant_count=$(bcftools view -H "$germline_vcf" | grep -c "|${gene}|")
                fi
                gene_counts+=("$variant_count")
                patient_total=$((patient_total + variant_count))
            done
            
            # Only show patients with variants
            if [ "$patient_total" -gt 0 ]; then
                printf "%-12s" "$patient_id" >> "$summary_file"
                for count in "${gene_counts[@]}"; do
                    printf "%-8s" "$count" >> "$summary_file"
                done
                printf "%-8s\n" "$patient_total" >> "$summary_file"
            fi
        fi
    fi
done

# Create high-impact variant summary
echo "" >> "$summary_file"
echo "=== HIGH-IMPACT VARIANTS BY TYPE ===" >> "$summary_file"
printf "%-25s %-8s %-20s\n" "Consequence" "Count" "Genes_Affected" >> "$summary_file"
echo "=================================================" >> "$summary_file"

# Count consequences from our collected data
for patient_dir in "$output_dir"/*/; do
    if [ -d "$patient_dir" ]; then
        patient_id=$(basename "$patient_dir")
        germline_vcf="${patient_dir}/final_consensus_pathogenic_germline_${patient_id}.vcf.gz"
        
        if [ -f "$germline_vcf" ]; then
            for gene in "${cancer_genes[@]}"; do
                bcftools view -H "$germline_vcf" | grep "|${gene}|" | while read line; do
                    info=$(echo "$line" | cut -f8)
                    
                    consequence="unknown"
                    
                    # Extract consequence from ANN or CSQ
                    if [[ "$info" =~ ANN=([^;]+) ]]; then
                        ann="${BASH_REMATCH[1]}"
                        IFS=',' read -ra ann_parts <<< "$ann"
                        for ann_part in "${ann_parts[@]}"; do
                            if [[ "$ann_part" =~ \|${gene}\| ]]; then
                                IFS='|' read -ra fields <<< "$ann_part"
                                consequence="${fields[1]:-unknown}"
                                break
                            fi
                        done
                    elif [[ "$info" =~ CSQ=([^;]+) ]]; then
                        csq="${BASH_REMATCH[1]}"
                        IFS=',' read -ra csq_parts <<< "$csq"
                        for csq_part in "${csq_parts[@]}"; do
                            if [[ "$csq_part" =~ \|${gene}\| ]]; then
                                IFS='|' read -ra fields <<< "$csq_part"
                                consequence="${fields[1]:-unknown}"
                                break
                            fi
                        done
                    fi
                    
                    # Count high-impact consequences
                    if [[ "$consequence" =~ stop_gained|frameshift|splice_acceptor|splice_donor|start_lost|stop_lost ]]; then
                        echo "$consequence:$gene" >> "/tmp/consequence_temp_$$"
                    fi
                done
            done
        fi
    fi
done

# Process consequence counts
if [ -f "/tmp/consequence_temp_$$" ]; then
    # Count occurrences and collect unique genes per consequence
    declare -A cons_counts
    declare -A cons_genes
    
    while read line; do
        consequence=$(echo "$line" | cut -d':' -f1)
        gene=$(echo "$line" | cut -d':' -f2)
        
        cons_counts["$consequence"]=$((${cons_counts["$consequence"]:-0} + 1))
        
        if [[ ! "${cons_genes["$consequence"]}" =~ $gene ]]; then
            if [ -z "${cons_genes["$consequence"]}" ]; then
                cons_genes["$consequence"]="$gene"
            else
                cons_genes["$consequence"]="${cons_genes["$consequence"]},$gene"
            fi
        fi
    done < "/tmp/consequence_temp_$$"
    
    # Output consequence summary
    for consequence in stop_gained frameshift_variant splice_acceptor_variant splice_donor_variant start_lost stop_lost; do
        count=${cons_counts[$consequence]:-0}
        genes="${cons_genes[$consequence]:-none}"
        if [ "$count" -gt 0 ]; then
            # Remove duplicates from gene list
            unique_genes=$(echo "$genes" | tr ',' '\n' | sort | uniq | tr '\n' ',' | sed 's/,$//')
            printf "%-25s %-8s %-20s\n" "$consequence" "$count" "$unique_genes" >> "$summary_file"
        fi
    done
    
    rm -f "/tmp/consequence_temp_$$"
else
    echo "No high-impact consequences found" >> "$summary_file"
fi

echo ""
echo "Analysis complete!"
echo ""
echo "Files created:"
echo "  Summary: $summary_file"
echo "  Detailed CSV: $detailed_file"
echo ""
echo "Key findings:"
echo "============="

# Display key findings to console
echo "Total patients analyzed: $total_patients_analyzed"
echo "Patients with cancer gene variants: $total_patients_with_variants"
if [ "$total_patients_analyzed" -gt 0 ]; then
    percentage=$(( (total_patients_with_variants * 100) / total_patients_analyzed ))
    echo "Percentage with cancer gene variants: $percentage%"
fi
echo ""

echo "Variants by gene:"
for gene in "${cancer_genes[@]}"; do
    variant_count=${gene_variant_counts[$gene]:-0}
    patient_list="${gene_patient_counts[$gene]:-}"
    patient_count=$(echo "$patient_list" | wc -w)
    
    if [ "$patient_count" -gt 0 ]; then
        echo "  $gene: $variant_count variants in $patient_count patients"
    else
        echo "  $gene: No variants found"
    fi
done

echo ""
echo "Frequency data status:"
for status in "no_population_data" "rare_variant" "common_variant"; do
    count=${frequency_status_counts[$status]:-0}
    if [ "$count" -gt 0 ]; then
        echo "  $status: $count variants"
    fi
done

echo ""
echo "Note: Many pathogenic cancer variants lack population frequency data"
echo "This is expected for rare/novel pathogenic mutations"
echo ""
echo "Check the summary file for detailed patient-by-gene breakdown!"
echo "Check the CSV file for variant-level details!"