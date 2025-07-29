#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process VCF2MAF_SOMATIC_WORKING {
    publishDir "${params.outdir}/somatic_mafs", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf_gz)
    
    output:
    tuple val(meta), path("${meta.id}_somatic.maf"), emit: maf
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}_somatic"
    
    """
    #!/bin/bash
    set +e
    set +u
    export BASHRCSOURCED=1
    
    # Set up environment
    source ~/.bashrc
    export PATH=/hpc/packages/minerva-rocky9/anaconda3/2024.06/bin:\$PATH
    source /hpc/packages/minerva-rocky9/anaconda3/2024.06/etc/profile.d/conda.sh
    conda activate vcf2maf_work
    module load samtools
    
    set -e
    set -u
    
    # Decompress VCF
    zcat ${vcf_gz} > input.vcf
    
    # Run vcf2maf
    /hpc/users/monsok03/software/vcf2maf/vcf2maf_fixed.pl \\
        --input-vcf input.vcf \\
        --output-maf ${prefix}.maf \\
        --tumor-id ${meta.id} \\
        --vcf-tumor-id ${meta.id} \\
        --species homo_sapiens \\
        --ncbi-build GRCh38 \\
        --vep-path /hpc/users/monsok03/.conda/envs/vcf2maf_work/bin \\
        --vep-data /sc/arion/projects/NGSCRC/Resources/VEP \\
        --cache-version 113 \\
        --ref-fasta /sc/arion/projects/NGSCRC/Resources/gatk_hg38/gatk_bundle_2024/Homo_sapiens_assembly38.fasta \\
        --samtools-exec \$(which samtools) \\
        --retain-info CSQ \\
        --vep-forks 1
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: 1.6.21
    END_VERSIONS
    """
}

process VCF2MAF_GERMLINE_WORKING {
    publishDir "${params.outdir}/germline_mafs", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf_gz)
    
    output:
    tuple val(meta), path("${meta.id}_germline.maf"), emit: maf
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}_germline"
    
    """
    #!/bin/bash
    set +e
    set +u
    export BASHRCSOURCED=1
    
    # Set up environment
    source ~/.bashrc
    export PATH=/hpc/packages/minerva-rocky9/anaconda3/2024.06/bin:\$PATH
    source /hpc/packages/minerva-rocky9/anaconda3/2024.06/etc/profile.d/conda.sh
    conda activate vcf2maf_work
    module load samtools
    
    set -e
    set -u
    
    # Decompress VCF
    zcat ${vcf_gz} > input.vcf
    
    # Run vcf2maf
    /hpc/users/monsok03/software/vcf2maf/vcf2maf_fixed.pl \\
        --input-vcf input.vcf \\
        --output-maf ${prefix}.maf \\
        --tumor-id ${meta.id} \\
        --normal-id ${meta.id} \\
        --species homo_sapiens \\
        --ncbi-build GRCh38 \\
        --vep-path /hpc/users/monsok03/.conda/envs/vcf2maf_work/bin \\
        --vep-data /sc/arion/projects/NGSCRC/Resources/VEP \\
        --cache-version 113 \\
        --ref-fasta /sc/arion/projects/NGSCRC/Resources/gatk_hg38/gatk_bundle_2024/Homo_sapiens_assembly38.fasta \\
        --samtools-exec \$(which samtools) \\
        --retain-info CSQ \\
        --vep-forks 1
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: 1.6.21
    END_VERSIONS
    """
}

process MERGE_ALL_MAFS {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path somatic_mafs
    path germline_mafs
    
    output:
    path "merged_consensus_all_samples.maf", emit: merged_maf
    path "merge_summary.txt", emit: summary
    
    script:
    """
    #!/bin/bash
    
    echo "=== Merging All MAF Files ===" > merge_summary.txt
    echo "Date: \$(date)" >> merge_summary.txt
    
    # Count input files
    SOMATIC_COUNT=\$(ls -1 ${somatic_mafs} 2>/dev/null | wc -l)
    GERMLINE_COUNT=\$(ls -1 ${germline_mafs} 2>/dev/null | wc -l)
    
    echo "Somatic MAF files: \$SOMATIC_COUNT" >> merge_summary.txt
    echo "Germline MAF files: \$GERMLINE_COUNT" >> merge_summary.txt
    
    # List all input files
    echo "" >> merge_summary.txt
    echo "Somatic MAF files:" >> merge_summary.txt
    ls -la ${somatic_mafs} >> merge_summary.txt 2>/dev/null || echo "No somatic MAFs" >> merge_summary.txt
    
    echo "" >> merge_summary.txt
    echo "Germline MAF files:" >> merge_summary.txt
    ls -la ${germline_mafs} >> merge_summary.txt 2>/dev/null || echo "No germline MAFs" >> merge_summary.txt
    
    # Create merged MAF file
    echo "Creating merged MAF file..." >> merge_summary.txt
    
    # Start with header from first file
    FIRST_MAF=\$(ls ${somatic_mafs} ${germline_mafs} 2>/dev/null | head -1)
    if [ -n "\$FIRST_MAF" ]; then
        head -2 "\$FIRST_MAF" > merged_consensus_all_samples.maf
        
        # Add all somatic variants (skip headers)
        for maf in ${somatic_mafs}; do
            if [ -f "\$maf" ]; then
                echo "Adding somatic variants from: \$maf" >> merge_summary.txt
                tail -n +3 "\$maf" >> merged_consensus_all_samples.maf
            fi
        done
        
        # Add all germline variants (skip headers)
        for maf in ${germline_mafs}; do
            if [ -f "\$maf" ]; then
                echo "Adding germline variants from: \$maf" >> merge_summary.txt
                tail -n +3 "\$maf" >> merged_consensus_all_samples.maf
            fi
        done
        
        # Count final variants
        TOTAL_VARIANTS=\$(tail -n +3 merged_consensus_all_samples.maf | wc -l)
        echo "" >> merge_summary.txt
        echo "Total variants in merged file: \$TOTAL_VARIANTS" >> merge_summary.txt
        echo "Merged MAF file size: \$(stat -c%s merged_consensus_all_samples.maf) bytes" >> merge_summary.txt
        
    else
        echo "ERROR: No MAF files found to merge!" >> merge_summary.txt
        touch merged_consensus_all_samples.maf
    fi
    
    echo "Merge completed: \$(date)" >> merge_summary.txt
    """
}

workflow {
    // Read and process the samplesheet
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id]
            def somatic_vcf = file(row.somatic_vcf, checkIfExists: true)
            def germline_vcf = (row.germline_vcf && row.germline_vcf != '' && row.germline_vcf != 'null') ? 
                               file(row.germline_vcf, checkIfExists: true) : null
            
            return [meta, somatic_vcf, germline_vcf]
        }
        .set { samples_ch }
    
    // Process somatic VCFs
    somatic_ch = samples_ch.map { meta, somatic_vcf, germline_vcf -> [meta, somatic_vcf] }
    VCF2MAF_SOMATIC_WORKING(somatic_ch)
    
    // Process germline VCFs (only for samples that have them)
    germline_ch = samples_ch
        .filter { meta, somatic_vcf, germline_vcf -> germline_vcf != null }
        .map { meta, somatic_vcf, germline_vcf -> [meta, germline_vcf] }
    
    VCF2MAF_GERMLINE_WORKING(germline_ch)
    
    // Merge all MAF files
    MERGE_ALL_MAFS(
        VCF2MAF_SOMATIC_WORKING.out.maf.map { meta, maf -> maf }.collect(),
        VCF2MAF_GERMLINE_WORKING.out.maf.map { meta, maf -> maf }.collect()
    )
}