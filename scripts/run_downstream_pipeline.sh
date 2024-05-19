#!/bin/bash

# Display message indicating post-upstream processing is starting
echo "Starting post-upstream processing..."

# Run post-upstream processing Python script
python scripts/post-upstream-processing.py data/metadata.csv data/counts_data.txt data/count_data.csv

# Display message indicating post-upstream processing is done
echo "Post-upstream processing done."

# Display message indicating differential gene expression analysis is starting
echo "Starting differential gene expression analysis..."

# Run differential gene expression analysis script
Rscript scripts/differential-gene-expression-analysis.R --input_count_data_path "data/count_data.csv" --metadata_path "data/metadata.csv"

# Display message indicating differential gene expression analysis is done
echo "Differential gene expression analysis done."

# Display message indicating gene set enrichment analysis is starting
echo "Starting gene set enrichment analysis..."

# Run gene set enrichment analysis script
Rscript scripts/gene-set-enrichment-analysis.R --deseq2_results_path data/dge/Deseq2-results-all.tsv --geneNamingTYpe "ENSEMBL"

# Display message indicating gene set enrichment analysis is done
echo "Gene set enrichment analysis done."



