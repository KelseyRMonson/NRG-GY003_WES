#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { VCF2MAF_SOMATIC } from './modules/vcf2maf_somatic'
include { VCF2MAF_GERMLINE } from './modules/vcf2maf_germline'
include { MERGE_MAFS } from './modules/merge_mafs'

workflow {
    // Read and process the samplesheet
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            println "Processing sample: ${row.sample_id}"
            println "  Somatic VCF: ${row.somatic_vcf}"
            println "  Germline VCF: ${row.germline_vcf ?: 'None'}"

            def meta = [id: row.sample_id]
            def somatic_vcf = file(row.somatic_vcf, checkIfExists: true)
            def germline_vcf = (row.germline_vcf && row.germline_vcf != '' && row.germline_vcf != 'null') ?
                               file(row.germline_vcf, checkIfExists: true) : null

            println "  Has germline: ${germline_vcf != null}"

            return [meta, somatic_vcf, germline_vcf]
        }
        .set { samples_ch }

    // Process somatic VCFs
    somatic_ch = samples_ch.map { meta, somatic_vcf, germline_vcf -> [meta, somatic_vcf] }
    VCF2MAF_SOMATIC(somatic_ch)  // Remove the reference_fasta parameter

    // Process germline VCFs (only for samples that have them)
    germline_ch = samples_ch
        .filter { meta, somatic_vcf, germline_vcf -> germline_vcf != null }
        .map { meta, somatic_vcf, germline_vcf -> [meta, germline_vcf] }

    VCF2MAF_GERMLINE(germline_ch)  // Remove the reference_fasta parameter

    // Merge MAFs
    MERGE_MAFS(
        VCF2MAF_SOMATIC.out.maf.collect(),
        VCF2MAF_GERMLINE.out.maf.collect()
    )
}