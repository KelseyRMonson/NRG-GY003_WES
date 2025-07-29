process MERGE_MAFS {
    publishDir "${params.outdir}/consensus", mode: 'copy'
    
    input:
    path somatic_mafs
    path germline_mafs
    
    output:
    path "all_samples_consensus.maf", emit: consensus_maf
    path "summary_stats.txt", emit: summary
    
    script:
    """
    # Create temporary files with mutation status annotations
    echo "Processing somatic MAFs..."
    SOMATIC_COUNT=0
    for maf in ${somatic_mafs}; do
        if [[ -s "\$maf" ]]; then
            awk 'BEGIN{FS=OFS="\\t"} NR==1 {print \$0, "Mutation_Status"} NR>1 {print \$0, "Somatic"}' "\$maf" > "\${maf%.maf}_with_status.maf"
            SOMATIC_COUNT=\$((SOMATIC_COUNT + 1))
        fi
    done
    
    echo "Processing germline MAFs..."
    GERMLINE_COUNT=0
    if [[ "${germline_mafs}" != "" ]]; then
        for maf in ${germline_mafs}; do
            if [[ -s "\$maf" ]]; then
                awk 'BEGIN{FS=OFS="\\t"} NR==1 {print \$0, "Mutation_Status"} NR>1 {print \$0, "Germline"}' "\$maf" > "\${maf%.maf}_with_status.maf"
                GERMLINE_COUNT=\$((GERMLINE_COUNT + 1))
            fi
        done
    fi
    
    # Find first non-empty file for header
    HEADER_FILE=""
    for file in *_with_status.maf; do
        if [[ -s "\$file" ]]; then
            HEADER_FILE="\$file"
            break
        fi
    done
    
    if [[ -z "\$HEADER_FILE" ]]; then
        echo "Error: No valid MAF files found!"
        exit 1
    fi
    
    # Start with header
    head -n1 "\$HEADER_FILE" > all_samples_consensus.maf
    
    # Add all data lines
    for file in *_with_status.maf; do
        if [[ -s "\$file" ]]; then
            tail -n +2 "\$file" >> all_samples_consensus.maf
        fi
    done
    
    # Generate summary statistics
    TOTAL_MUTATIONS=\$(tail -n +2 all_samples_consensus.maf | wc -l)
    SOMATIC_MUTATIONS=\$(tail -n +2 all_samples_consensus.maf | awk -F'\\t' '\$NF=="Somatic"' | wc -l)
    GERMLINE_MUTATIONS=\$(tail -n +2 all_samples_consensus.maf | awk -F'\\t' '\$NF=="Germline"' | wc -l)
    UNIQUE_SAMPLES=\$(tail -n +2 all_samples_consensus.maf | cut -f16 | sort -u | wc -l)
    
    cat > summary_stats.txt <<EOF
VCF to MAF Conversion Summary
=============================
Total mutations: \$TOTAL_MUTATIONS
Somatic mutations: \$SOMATIC_MUTATIONS
Germline mutations: \$GERMLINE_MUTATIONS
Unique samples: \$UNIQUE_SAMPLES
Somatic MAF files processed: \$SOMATIC_COUNT
Germline MAF files processed: \$GERMLINE_COUNT
EOF
    
    echo "Summary:"
    cat summary_stats.txt
    """
    
    stub:
    """
    touch all_samples_consensus.maf
    touch summary_stats.txt
    """
}