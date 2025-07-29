#!/bin/bash

# Clean up previous run
rm -rf work/
rm -rf results/
rm -rf .nextflow*

# Clean up trace files
mkdir -p results/pipeline_info

echo "Cleaned up previous run. Ready to restart."