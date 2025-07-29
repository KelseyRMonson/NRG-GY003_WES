# NRG-GY003_WES
Whole Exome Sequencing analysis from NRG-GY003 clinical trial

- [01_Sarek](01_Sarek) Folder
    - Folder containing initial Nextflow [Sarek](https://nf-co.re/sarek/3.4.0/) pipeline output
    - The Sarek pipeline was used to process raw .fastq files into annotated vcfs
    - This folder contains
        - MultiQC report and plots
        - Pipeline info
        - *Does not contain actual variant called files for size and privacy reasons*
- [02_vcf2maf](02_vcf2maf) Folder
    - Folder containing post-processing of vcfs from Sarek to consensus maf file for downstream analysis
    -  Generated a consensus maf using the [vcf2maf](https://nf-co.re/modules/vcf2maf) Nextflow pipeline for all samples, annotating variants as either somatic or germline.
<details>
    <summary>Variant processing details</summary>
    
- Variants were processed as follows:
    - Find consensus calls between >= 2 variant callers
    - Drop variants without PASS status
    - Further refine PASS variants to those with minimum tumor depth of 20, and an allele frequency 0.05 < AF < 0.95
    - Identify rare pathogenic germline variants
        -  Protein-coding variants annotated with HIGH impact
        -  Present in >= 2 callers
        -  <5% population frequency (likely a bit liberal; "rare" variants typically are <1%, may refine further if needed)
        -  Include all likely oncogenes regardless of population frequency ("BRCA1" "BRCA2" "TP53" "PTEN" "ATM" "CHEK2" "PALB2" "APC" "MUTYH")
-  This generated two vcfs, one with the processed somatic variants passing filters and present in >=2 callers, and one with likely pathogenic germline variants

</details>

- [03_maftools](03_maftools) Folder 🚧 
    - Folder containing R Project files for maftools analysis
    - *Work in progress*
