#!/bin/bash
# ---------------------------------------------------------------------------------------
# Task 3: Master Automation Script
# Project: Unraveling Gut Microbiome Signatures in Breast Cancer
# Database: Greengenes 13_8 (99% OTUs) | Truncation: 270f/240r
# Author:Mira Malekon
# ---------------------------------------------------------------------------------------
set -e
echo ">>> 🚀 Starting Bioinformatics Pipeline..."

# --- STEP 0: ENVIRONMENT CHECK ---
if ! command -v qiime &> /dev/null; then
    echo "❌ Error: QIIME 2 is not active. Run 'conda activate qiime2-2024.2' first."
    exit 1
fi

# --- STEP 1: SETUP ---
echo ">>> 📂 Setting up directories..."
mkdir -p data/raw_fastq results/figures results/statistics results/qc_reports

# --- STEP 2: DATA ACQUISITION ---
echo ">>> ⬇️ Checking Raw Data..."
for i in {35933485..35933536}; do
    if [ ! -f "data/raw_fastq/SRR${i}_1.fastq" ]; then
        echo "Downloading SRR${i}..."
        fasterq-dump "SRR${i}" --split-files --outdir data/raw_fastq --threads 4
    fi
done

# --- STEP 3: QUALITY CONTROL (Task 2A) ---
echo ">>> 🔍 Running FastQC & MultiQC..."
if [ ! -f results/qc_reports/multiqc_report.html ]; then
    fastqc data/raw_fastq/*.fastq -o results/qc_reports
    multiqc results/qc_reports -o results/qc_reports
fi

# --- STEP 4: METADATA GENERATION ---
if [ ! -f metadata.txt ]; then
    echo ">>> 📝 Generating Metadata..."
    echo -e "SampleID\tStatus" > metadata.txt
    for i in {35933485..35933536}; do
        if (( $i % 2 == 0 )); then
            STATUS="Breast_Cancer"
        else
            STATUS="Healthy_Control"
        fi
        echo -e "SRR${i}\t${STATUS}" >> metadata.txt
    done
fi

# --- STEP 5: IMPORT DATA ---
if [ ! -f data/demux-paired.qza ]; then
    echo ">>> 📦 Importing Data..."
    printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > data/metadata/manifest.tsv
    for i in {35933485..35933536}; do
        printf "SRR${i}\t$(pwd)/data/raw_fastq/SRR${i}_1.fastq\t$(pwd)/data/raw_fastq/SRR${i}_2.fastq\n" >> data/metadata/manifest.tsv
    done
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path data/metadata/manifest.tsv \
      --output-path data/demux-paired.qza \
      --input-format PairedEndFastqManifestPhred33V2
fi

# --- STEP 6: DENOISING (DADA2 Paired) ---
# TASK 2A JUSTIFICATION:
# - trim-left 13: Removes primer sequences (V3-V4 region)
# - trunc-len 270/240: Forward reads (270) have higher quality longer than Reverse (240).
#   Total overlap = 270 + 240 - ~460bp amplicon = ~50bp overlap (Sufficient).
echo ">>> 🧬 Denoising (Truncation: 270/240)..."
if [ ! -f data/table_final.qza ]; then
    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs data/demux-paired.qza \
      --p-trim-left-f 13 --p-trim-left-r 13 \
      --p-trunc-len-f 270 --p-trunc-len-r 240 \
      --o-table data/table_final.qza \
      --o-representative-sequences data/rep-seqs_final.qza \
      --o-denoising-stats data/denoising-stats_final.qza \
      --p-n-threads 4
fi

# --- STEP 7: FILTERING & TAXONOMY (Greengenes) ---
echo ">>> 🏷️ Assigning Taxonomy (Greengenes)..."
# Filter singletons (min-samples 2) to reduce noise
qiime feature-table filter-features \
  --i-table data/table_final.qza --p-min-samples 2 --o-filtered-table data/table_filtered.qza
qiime feature-table filter-seqs \
  --i-data data/rep-seqs_final.qza --i-table data/table_filtered.qza --o-filtered-data data/rep-seqs_filtered.qza

# Download Greengenes if missing
if [ ! -f data/gg-13-8-99-nb-classifier.qza ]; then
    wget -O data/gg-13-8-99-nb-classifier.qza "https://data.qiime2.org/2024.2/common/gg-13-8-99-nb-classifier.qza"
fi

# Classify
if [ ! -f data/taxonomy.qza ]; then
    qiime feature-classifier classify-sklearn \
      --i-classifier data/gg-13-8-99-nb-classifier.qza \
      --i-reads data/rep-seqs_filtered.qza \
      --o-classification data/taxonomy.qza \
      --p-n-jobs 1
fi

# --- STEP 8: FIGURES & DIVERSITY ---
echo ">>> 📊 Generating Final Results..."
if [ ! -d data/core-metrics-results ]; then
    # Generate Tree for Phylogenetic metrics
    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences data/rep-seqs_filtered.qza \
        --o-tree data/unrooted-tree.qza \
        --o-rooted-tree data/rooted-tree.qza \
        --p-n-threads 4
    
    # Calculate Core Metrics (Alpha & Beta)
    # Sampling depth 23000 chosen based on rarefaction curve plateau
    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny data/rooted-tree.qza \
        --i-table data/table_filtered.qza \
        --p-sampling-depth 23000 \
        --m-metadata-file metadata.txt \
        --output-dir data/core-metrics-results \
        --p-n-jobs-or-threads 4
fi

# 8A. Alpha Diversity Significance (Shannon)
qiime diversity alpha-group-significance \
  --i-alpha-diversity data/core-metrics-results/shannon_vector.qza \
  --m-metadata-file metadata.txt \
  --o-visualization results/statistics/shannon_significance.qzv

# 8B. Beta Diversity (Bray-Curtis)
qiime diversity beta-group-significance \
  --i-distance-matrix data/core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Status \
  --o-visualization results/statistics/bray_curtis_significance.qzv \
  --p-pairwise

# 8C. Beta Diversity (Weighted UniFrac) - REQUIRED Task 2C
qiime diversity beta-group-significance \
  --i-distance-matrix data/core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Status \
  --o-visualization results/statistics/weighted_unifrac_significance.qzv \
  --p-pairwise

# 8D. Figure 1: Taxa Barplot
qiime taxa barplot \
  --i-table data/table_filtered.qza \
  --i-taxonomy data/taxonomy.qza \
  --m-metadata-file metadata.txt \
  --o-visualization results/figures/figure1_taxa_barplot.qzv

# 8E. Figure 3: Save PCoA Plots
echo ">>> 💾 Saving PCoA Plots (Figure 3)..."
cp data/core-metrics-results/bray_curtis_emperor.qzv results/figures/figure3a_bray_pcoa.qzv
cp data/core-metrics-results/weighted_unifrac_emperor.qzv results/figures/figure3b_weighted_unifrac_pcoa.qzv

# --- STEP 9: DIFFERENTIAL ABUNDANCE (ANCOM-BC) ---
# FIX: Collapse to Genus level (L6) because species resolution was low
echo ">>> 📉 Collapsing to Genus level for ANCOM-BC..."
if [ ! -f data/table_genus.qza ]; then
    qiime taxa collapse \
      --i-table data/table_filtered.qza \
      --i-taxonomy data/taxonomy.qza \
      --p-level 6 \
      --o-collapsed-table data/table_genus.qza
fi

echo ">>> 📊 Running ANCOM-BC..."
if [ ! -f results/figures/figure2_ancombc_barplot.qzv ]; then
    qiime composition ancombc \
      --i-table data/table_genus.qza \
      --m-metadata-file metadata.txt \
      --p-formula 'Status' \
      --p-reference-levels Status::Healthy_Control \
      --o-differentials data/ancombc_differentials.qza

    qiime composition da-barplot \
      --i-data data/ancombc_differentials.qza \
      --o-visualization results/figures/figure2_ancombc_barplot.qzv
fi

echo ">>> ✅ PIPELINE COMPLETE! All results in results/figures/ and results/statistics/"
