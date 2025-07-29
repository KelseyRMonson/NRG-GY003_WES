process VCF2MAF_SOMATIC {
    publishDir "${params.outdir}/somatic_mafs", mode: 'copy'

    input:
    tuple val(meta), path(vcf_gz)

    output:
    tuple val(meta), path("${meta.id}_somatic.maf"), emit: maf
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_somatic"

    """
    #!/bin/bash
    set +e
    set +u
    export BASHRCSOURCED=1

    # Set up environment (exactly like the working test)
    source ~/.bashrc
    export PATH=/hpc/packages/minerva-rocky9/anaconda3/2024.06/bin:\$PATH
    source /hpc/packages/minerva-rocky9/anaconda3/2024.06/etc/profile.d/conda.sh
    conda activate vcf2maf_work
    module load samtools

    set -e
    set -u

    # Decompress VCF
    zcat ${vcf_gz} > input.vcf

    # Run vcf2maf (exactly like the working test)
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
        vcf2maf: \$(vcf2maf_fixed.pl --help | grep -i version | head -1 | sed 's/.*version //i' || echo "1.6.21")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_somatic"
    """
    touch ${prefix}.maf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: 1.6.21
    END_VERSIONS
    """
}
