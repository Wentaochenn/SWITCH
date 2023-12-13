# SWITCH: a robust and selective whole-genome amplification tool for Chlamydia trachomatis

SWITCH is a excellent tool designed for robust and selective whole-genome amplification, advancing molecular surveillance of Chlamydia trachomatis.

## Detailed Tutorials for SWITCH Scripts
### 1. reads_processing_pipeline.py

#### Purpose
This script is a pipeline that generates pseudosequences from raw sequencing data. The pipeline is designed to be run using Snakemake, a workflow management system that automates the execution of complex data processing.

### 2. Genome_depth_plotter.py
#### Purpose
`Genome_depth_plotter.py` is used to visualize the genome coverage of Chlamydia trachomatis from sequencing data. The input is a BAM file, and the output is a graphical representation of the genome coverage.

### 3. Primer_match_calculator.py

#### Purpose
This script calculates the targeting efficiency of primer sets designed for SWITCH. It performs in-silico analysis to assess how effectively the primers target the Chlamydia trachomatis and human genomes.

