#!/bin/bash
#BSUB -J vcf2maf_nextflow_head
#BSUB -P acc_NGSCRC
#BSUB -q premium
#BSUB -n 1
#BSUB -W 12:00
#BSUB -R rusage[mem=8000]
#BSUB -o nextflow_head_%J.out
#BSUB -e nextflow_head_%J.err
#BSUB -L /bin/bash

# Exit on any error
set -euo pipefail

# Set proxy settings
export http_proxy=http://172.28.7.1:3128
export https_proxy=http://172.28.7.1:3128
export all_proxy=http://172.28.7.1:3128

# Reduce Java memory usage for the head job
export JAVA_OPTS="-Xms1g -Xmx4g -XX:+UseSerialGC"
export NXF_OPTS="-Xms1g -Xmx4g -XX:+UseSerialGC"

# Disable virtual threads to reduce threading issues
export NXF_ENABLE_VIRTUAL_THREADS=false

# Set working directory
cd /sc/arion/projects/NGSCRC/Scripts/NRG-GY003/WES/manage_Sarek_output/vcf2maf_pipeline

# Load required modules
module load java/21.0.4
module load nextflow

echo "=== Starting VCF2MAF Pipeline ==="
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo "Java version: $(java -version 2>&1 | head -1)"
echo "Nextflow version: $(nextflow -version)"

# Clean up previous run
echo "Cleaning up previous runs..."
rm -rf work/ .nextflow* final_results/

# Create samplesheet
echo "Creating samplesheet..."
if [ ! -f "create_samplesheet.sh" ]; then
    echo "ERROR: create_samplesheet.sh not found!"
    exit 1
fi

chmod +x create_samplesheet.sh
./create_samplesheet.sh

# Verify samplesheet was created
if [ ! -f "samplesheet.csv" ]; then
    echo "ERROR: samplesheet.csv was not created!"
    exit 1
fi

echo "Samplesheet created with $(tail -n +2 samplesheet.csv | wc -l) samples"
echo "First few lines of samplesheet:"
head -5 samplesheet.csv

# Verify pipeline file exists
if [ ! -f "complete_working_pipeline.nf" ]; then
    echo "ERROR: complete_working_pipeline.nf not found!"
    exit 1
fi

# Create output directory
mkdir -p final_results/pipeline_info

echo "=== Starting Nextflow Pipeline ==="
echo "Date: $(date)"

# Run the pipeline with error handling
nextflow run complete_working_pipeline.nf \
    -profile lsf \
    --input samplesheet.csv \
    --outdir final_results \
    -with-report final_results/pipeline_info/execution_report.html \
    -with-timeline final_results/pipeline_info/execution_timeline.html \
    -with-trace final_results/pipeline_info/execution_trace.txt \
    -with-dag final_results/pipeline_info/pipeline_dag.svg \
    -resume

PIPELINE_EXIT_CODE=$?

echo "=== Pipeline Completed ==="
echo "Date: $(date)"
echo "Exit code: $PIPELINE_EXIT_CODE"

# Check results
if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo "=== Pipeline SUCCESS! ==="
    echo "Results summary:"
    echo "Somatic MAFs: $(ls final_results/somatic_mafs/*.maf 2>/dev/null | wc -l)"
    echo "Germline MAFs: $(ls final_results/germline_mafs/*.maf 2>/dev/null | wc -l)"
    echo "Merged MAF: $(ls final_results/merged_consensus_all_samples.maf 2>/dev/null | wc -l)"
    
    if [ -f "final_results/merged_consensus_all_samples.maf" ]; then
        echo "Merged MAF variants: $(tail -n +3 final_results/merged_consensus_all_samples.maf | wc -l)"
    fi
else
    echo "=== Pipeline FAILED! ==="
    echo "Check the Nextflow log and work directories for details"
    echo "Recent work directories:"
    find work/ -name ".command.log" -newer work 2>/dev/null | head -5 | xargs dirname || echo "No work directories found"
fi

exit $PIPELINE_EXIT_CODE