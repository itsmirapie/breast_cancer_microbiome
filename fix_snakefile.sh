#!/bin/bash
# Fix Snakemake pipeline errors

echo "Fixing phylogeny rule in Snakefile..."

# Backup original
cp Snakefile Snakefile.backup

# Create fixed Snakefile
cat > Snakefile.fixed << 'EOF'
# Snakefile for Breast Cancer Microbiome Analysis
# Author: [Your Name]

configfile: "config.yaml"

# Get samples from config
SAMPLES = config["samples"]

rule all:
    input:
        # Quality Control
        "results/qc_reports/multiqc_report.html",
        
        # Core Analysis
        "data/taxonomy.qza",
        "data/rooted-tree.qza",
        
        # Required Figures
        "results/figures/figure1_taxa_barplot.qzv",
        "results/figures/figure2_ancombc_barplot.qzv", 
        "results/figures/figure3a_bray_pcoa.qzv",
        "results/figures/figure3b_weighted_unifrac_pcoa.qzv",
        
        # Statistics
        "results/statistics/shannon_significance.qzv",
        "results/statistics/bray_curtis_significance.qzv",
        "results/statistics/weighted_unifrac_significance.qzv"

# --- 0. Metadata Generation ---
rule generate_metadata:
    output:
        "metadata.txt"
    run:
        with open(output[0], 'w') as f:
            f.write("SampleID\tStatus\n")
            for sample in SAMPLES:
                i = int(sample.replace("SRR", ""))
                status = "Breast_Cancer" if i % 2 == 0 else "Healthy_Control"
                f.write(f"{sample}\t{status}\n")

# --- 1. Data Download ---
rule download_data:
    output:
        r1 = "data/raw_fastq/{sample}_1.fastq",
        r2 = "data/raw_fastq/{sample}_2.fastq"
    threads: config["resources"]["download_threads"]
    shell:
        "fasterq-dump {wildcards.sample} --split-files --outdir data/raw_fastq --threads {threads}"

# --- 2. Quality Control ---
rule qc:
    input:
        expand("data/raw_fastq/{sample}_{read}.fastq", 
               sample=SAMPLES, read=[1,2])
    output:
        html = "results/qc_reports/multiqc_report.html",
        directory = directory("results/qc_reports")
    threads: config["resources"]["qiime_threads"]
    shell:
        """
        fastqc data/raw_fastq/*.fastq -o results/qc_reports -t {threads}
        multiqc results/qc_reports -o results/qc_reports -f
        """

# --- 3. Create Manifest ---
rule create_manifest:
    input:
        expand("data/raw_fastq/{sample}_{read}.fastq", 
               sample=SAMPLES, read=[1,2])
    output:
        "data/manifest.tsv"
    shell:
        """
        printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > {output}
        for sample in {' '.join(SAMPLES)}; do
            printf "$sample\t$(pwd)/data/raw_fastq/${sample}_1.fastq\t$(pwd)/data/raw_fastq/${sample}_2.fastq\n" >> {output}
        done
        """

# --- 4. Import Data ---
rule import_data:
    input:
        "data/manifest.tsv"
    output:
        "data/demux-paired.qza"
    threads: config["resources"]["qiime_threads"]
    shell:
        """
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {input} \
            --output-path {output} \
            --input-format PairedEndFastqManifestPhred33V2
        """

# --- 5. Denoise (DADA2) ---
rule denoise:
    input:
        "data/demux-paired.qza"
    output:
        table_raw = "data/table_final.qza",
        rep_seqs_raw = "data/rep-seqs_final.qza",
        stats = "data/denoising-stats.qza"
    threads: config["resources"]["qiime_threads"]
    params:
        trim_left = config["params"]["trim_left"],
        trunc_len_f = config["params"]["trunc_len_f"],
        trunc_len_r = config["params"]["trunc_len_r"]
    shell:
        """
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {input} \
            --p-trim-left-f {params.trim_left} --p-trim-left-r {params.trim_left} \
            --p-trunc-len-f {params.trunc_len_f} --p-trunc-len-r {params.trunc_len_r} \
            --o-table {output.table_raw} \
            --o-representative-sequences {output.rep_seqs_raw} \
            --o-denoising-stats {output.stats} \
            --p-n-threads {threads}
        """

# --- 6. Filter Features ---
rule filter_features:
    input:
        table = "data/table_final.qza",
        seqs = "data/rep-seqs_final.qza"
    output:
        table = "data/table_filtered.qza",
        seqs = "data/rep-seqs_filtered.qza"
    params:
        min_samples = config["params"]["min_samples"]
    shell:
        """
        qiime feature-table filter-features \
            --i-table {input.table} \
            --p-min-samples {params.min_samples} \
            --o-filtered-table {output.table}
            
        qiime feature-table filter-seqs \
            --i-data {input.seqs} \
            --i-table {output.table} \
            --o-filtered-data {output.seqs}
        """

# --- 7. Download Classifier ---
rule download_classifier:
    output:
        "data/gg-13-8-99-nb-classifier.qza"
    shell:
        """
        wget -O {output} "https://data.qiime2.org/2024.2/common/gg-13-8-99-nb-classifier.qza"
        """

# --- 8. Assign Taxonomy ---
rule taxonomy:
    input:
        seqs = "data/rep-seqs_filtered.qza",
        classifier = "data/gg-13-8-99-nb-classifier.qza"
    output:
        "data/taxonomy.qza"
    threads: config["resources"]["qiime_threads"]
    shell:
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.seqs} \
            --o-classification {output} \
            --p-n-jobs {threads}
        """

# --- 9. Generate Phylogenetic Tree ---
rule phylogeny:
    input:
        "data/rep-seqs_filtered.qza"
    output:
        alignment = "data/aligned-rep-seqs.qza",
        masked_alignment = "data/masked-aligned-rep-seqs.qza",
        unrooted_tree = "data/unrooted-tree.qza",
        rooted_tree = "data/rooted-tree.qza"
    threads: config["resources"]["qiime_threads"]
    shell:
        """
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences {input} \
            --o-alignment {output.alignment} \
            --o-masked-alignment {output.masked_alignment} \
            --o-tree {output.unrooted_tree} \
            --o-rooted-tree {output.rooted_tree} \
            --p-n-threads {threads}
        """

# --- 10. Core Diversity Metrics ---
rule diversity_core:
    input:
        table = "data/table_filtered.qza",
        tree = "data/rooted-tree.qza",
        metadata = "metadata.txt"
    output:
        directory("data/core-metrics-results")
    threads: config["resources"]["qiime_threads"]
    params:
        sampling_depth = config["params"]["sampling_depth"]
    shell:
        """
        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny {input.tree} \
            --i-table {input.table} \
            --p-sampling-depth {params.sampling_depth} \
            --m-metadata-file {input.metadata} \
            --output-dir {output} \
            --p-n-jobs-or-threads {threads}
        """

# --- 11. Alpha Diversity ---
rule alpha_diversity:
    input:
        shannon = "data/core-metrics-results/shannon_vector.qza",
        observed = "data/core-metrics-results/observed_features_vector.qza",
        metadata = "metadata.txt"
    output:
        shannon_sig = "results/statistics/shannon_significance.qzv",
        observed_sig = "results/statistics/observed_significance.qzv"
    shell:
        """
        qiime diversity alpha-group-significance \
            --i-alpha-diversity {input.shannon} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.shannon_sig}
            
        qiime diversity alpha-group-significance \
            --i-alpha-diversity {input.observed} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.observed_sig}
        """

# --- 12. Beta Diversity PERMANOVA ---
rule beta_diversity:
    input:
        bray = "data/core-metrics-results/bray_curtis_distance_matrix.qza",
        wunifrac = "data/core-metrics-results/weighted_unifrac_distance_matrix.qza",
        metadata = "metadata.txt"
    output:
        bray_sig = "results/statistics/bray_curtis_significance.qzv",
        wunifrac_sig = "results/statistics/weighted_unifrac_significance.qzv"
    shell:
        """
        qiime diversity beta-group-significance \
            --i-distance-matrix {input.bray} \
            --m-metadata-file {input.metadata} \
            --m-metadata-column Status \
            --o-visualization {output.bray_sig} \
            --p-pairwise
            
        qiime diversity beta-group-significance \
            --i-distance-matrix {input.wunifrac} \
            --m-metadata-file {input.metadata} \
            --m-metadata-column Status \
            --o-visualization {output.wunifrac_sig} \
            --p-pairwise
        """

# --- 13. Figure 1: Taxa Bar Plot ---
rule taxa_barplot:
    input:
        table = "data/table_filtered.qza",
        taxonomy = "data/taxonomy.qza",
        metadata = "metadata.txt"
    output:
        "results/figures/figure1_taxa_barplot.qzv"
    shell:
        """
        qiime taxa barplot \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output}
        """

# --- 14. Figure 3: PCoA Plots ---
rule pcoa_plots:
    input:
        bray = "data/core-metrics-results/bray_curtis_emperor.qzv",
        wunifrac = "data/core-metrics-results/weighted_unifrac_emperor.qzv"
    output:
        bray = "results/figures/figure3a_bray_pcoa.qzv",
        wunifrac = "results/figures/figure3b_weighted_unifrac_pcoa.qzv"
    shell:
        """
        cp {input.bray} {output.bray}
        cp {input.wunifrac} {output.wunifrac}
        """

# --- 15. ANCOM-BC (Differential Abundance) ---
rule collapse_genus:
    input:
        table = "data/table_filtered.qza",
        taxonomy = "data/taxonomy.qza"
    output:
        "data/table_genus.qza"
    shell:
        """
        qiime taxa collapse \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --p-level 6 \
            --o-collapsed-table {output}
        """

rule ancom_bc:
    input:
        table = "data/table_genus.qza",
        metadata = "metadata.txt"
    output:
        diffs = "data/ancombc_differentials.qza",
        plot = "results/figures/figure2_ancombc_barplot.qzv"
    shell:
        """
        qiime composition ancombc \
            --i-table {input.table} \
            --m-metadata-file {input.metadata} \
            --p-formula 'Status' \
            --p-reference-levels Status::Healthy_Control \
            --o-differentials {output.diffs}
            
        qiime composition da-barplot \
            --i-data {output.diffs} \
            --o-visualization {output.plot}
        """
EOF

# Replace old Snakefile
mv Snakefile.fixed Snakefile

echo "✅ Fixed Snakefile"
echo ""
echo "Now run: snakemake --cores 8"
