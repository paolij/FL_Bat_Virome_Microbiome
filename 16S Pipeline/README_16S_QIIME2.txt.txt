# 16S rRNA Amplicon Analysis with QIIME 2

Reproducible pipeline for demultiplexed paired-end 16S reads, including primer trimming, read joining, quality filtering, Deblur denoising, taxonomic classification, and BIOM/TSV export.

Important:
- Create a separate SLURM job script for each step (import & cutadapt, deblur, tree/rarefaction, taxonomy, export).
- Keep all scripts and outputs in the same project folder for consistency and easier resumption.
- Use consistent artifact names (step.01a…, step.02.a…, etc.) across scripts.

---

## Prerequisites

- QIIME 2 available on your system (e.g., module qiime2/2024.10)
- Demultiplexed paired-end FASTQ files in Casava 1.8 format
- SILVA V3–V4 classifier (.qza)
- Sample metadata file (tab-delimited, UTF-8)

---

## Suggested Project Layout

- project/
  - input/ (demultiplexed FASTQs)
  - classifier/ (e.g., silva-v3v4-classifier.qza)
  - scripts/ (SLURM scripts, one per step)
  - outputs/ (all .qza/.qzv and exported files)

Place all SLURM scripts in project/scripts and set their working directory to project/outputs (or project/). Keep all outputs in one place.

---

## SLURM Usage (Template)

Create one SLURM script per step (e.g., 01_cutadapt.slurm, 02_deblur.slurm, 03_tree_rarefaction.slurm, 04_taxonomy.slurm, 05_export.slurm). Example template:

```bash
#!/bin/bash
#SBATCH --job-name=qiime_step
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

module load qiime2/2024.10

# Create and move to a consistent working directory for all steps
mkdir -p outputs logs
cd outputs

# User-configurable variables (keep generic across steps)
INPUT_DIR="path/to/input"
METADATA="path/to/metadata.tsv"
CLASSIFIER="path/to/silva_v3v4_classifier.qza"
FORWARD_PRIMER="CCTACGGGNGGCWGCAG"   # 314F
REVERSE_PRIMER="GACTACHVGGGTATCTAATCC" # 805R

# Step-specific commands go here...

## STEP 1 — Import demultiplexed paired-end reads

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --input-path "$INPUT_DIR" \
  --output-path step.01a.demux.pe.qza

# STEP 1b — Primer trimming (V3–V4: 314F/805R)
# 314F (forward): CCTACGGGNGGCWGCAG
# 805R (reverse): GACTACHVGGGTATCTAATCC
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences step.01a.demux.pe.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences step.01b.cutadapt.qza \
  --p-match-read-wildcards --verbose

# STEP 1c — Summarize post-trim demux
qiime demux summarize \
  --i-data step.01b.cutadapt.qza \
  --o-visualization step.01b.cutadapt.qzv

# NOTE: Join paired reads before quality filtering (artifact name expected below):
# Expected joined artifact: step.01c.joined.qza
# If you haven't joined reads yet, insert your join step here (e.g., q2-vsearch join-pairs).

# STEP 1d — Quality filtering (q-score)
# Set matplotlib backend to avoid display issues on headless systems
echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

qiime demux summarize \
  --i-data  step.01c.joined.qza \
  --o-visualization  step.01c.joined.qzv

qiime quality-filter q-score \
  --i-demux step.01c.joined.qza \
  --p-min-quality 4 \
  --p-quality-window 3 \
  --verbose \
  --o-filtered-sequences step.01d.joined.filtered.qza \
  --o-filter-stats step.01d.joined.filtered.stats.qza

qiime demux summarize \
  --i-data step.01d.joined.filtered.qza \
  --o-visualization step.01d.joined.filtered.qzv

date

## STEP 2 — Deblur denoise-16S
qiime deblur denoise-16S \
  --i-demultiplexed-seqs step.01d.joined.filtered.qza \
  --p-trim-length 400 \
  --p-mean-error 0.00125 \
  --p-min-reads 10 \
  --p-sample-stats \
  --o-representative-sequences step.02.a.deblur.repseqs.qza \
  --o-table step.02.a.deblur.table.qza \
  --o-stats step.02.a.deblur.stats.qza

# STEP 2b — Summaries
qiime feature-table summarize \
  --i-table step.02.a.deblur.table.qza \
  --m-sample-metadata-file "$METADATA" \
  --o-visualization step.02.a.deblur.summary.qzv

qiime feature-table tabulate-seqs \
  --i-data step.02.a.deblur.repseqs.qza \
  --o-visualization step.02.a.deblur.repseqs.qzv

date

## STEP 3 — Phylogeny and alpha rarefaction

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences step.02.a.deblur.repseqs.qza \
  --o-alignment step.03.a.deblur.aligned.repseqs.qza \
  --o-masked-alignment step.03.a.deblur.aligned.masked.repseqs.qza \
  --o-tree tree.unrooted.qza \
  --o-rooted-tree tree.rooted.qza

qiime diversity alpha-rarefaction \
  --i-table step.02.a.deblur.table.qza \
  --i-phylogeny tree.rooted.qza \
  --m-metadata-file "$METADATA" \
  --p-min-depth 1 \
  --p-max-depth 10000 \
  --p-metrics observed_features --p-metrics shannon --p-metrics chao1 --p-metrics faith_pd \
  --output-dir rarefactionCurves

date

## STEP 4 — Taxonomic classification

# Choosing sampling depth for core-metrics:

# Open rarefactionCurves/alpha-rarefaction.qzv in the QIIME 2 viewer and identify where curves approach a plateau.
# Select a depth that reaches the plateau while losing as few samples as possible; this will be your SAMP_DEPTH for core-metrics [1].

module load qiime2/2024.10

#Define variables of interest to be used later for the group comparisons of beta diversity:
sampDepth=10000                    		# set from rarefaction plateau
myCols=("group" "sample_type" "sample_id")      # update with your metadata columns

# Need to set the depth so you can reach plateu on the rarifaction curve, but lose as few samples as possible
# Run core-metrics-phylogenetic:
# Matrix used for either alpha or beta diversity will tell you the number of samples, number of sequences for each sample, as well as the genetic distance between each of the samples

# Step 4.1 Run core-metrics-phylogenetic:
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny tree.rooted.qza \
    --i-table step.02.a.deblur.table.qza \
    --p-sampling-depth $sampDepth \
    --m-metadata-file $metadata \
    --output-dir coreMetrics

# Step 4.2 Alpha diversity metrics significance
qiime diversity alpha-group-significance \
    --i-alpha-diversity coreMetrics/faith_pd_vector.qza \
    --m-metadata-file $metadata \
    --o-visualization coreMetrics/faith_pd.significance.qzv

qiime diversity alpha-group-significance \
    --i-alpha-diversity coreMetrics/shannon_vector.qza \
    --m-metadata-file $metadata \
    --o-visualization coreMetrics/shannon.significance.qzv


qiime diversity alpha-group-significance \
    --i-alpha-diversity coreMetrics/evenness_vector.qza \
    --m-metadata-file $metadata \
    --o-visualization coreMetrics/evenness.significance.qzv


# Step 4.3 Beta diversity metrics significance 
#Loop each variable:
nCols=${#myCols[@]}
for i in $(seq 0 $((nCols-1)))
do
    qiime diversity beta-group-significance \
        --i-distance-matrix coreMetrics/unweighted_unifrac_distance_matrix.qza \
        --m-metadata-file $metadata \
        --m-metadata-column ${myCols[$i]} \
        --o-visualization coreMetrics/unifrac.significance.${myCols[$i]}.qzv \
        --p-pairwise

    qiime diversity beta-group-significance \
        --i-distance-matrix coreMetrics/weighted_unifrac_distance_matrix.qza \
        --m-metadata-file $metadata \
        --m-metadata-column ${myCols[$i]} \
        --o-visualization coreMetrics/wunifrac.significance.${myCols[$i]}.qzv \
        --p-pairwise
done

# Step 4.4 Taxonomic classification
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads step.02.a.deblur.repseqs.qza \
  --o-classification step.04.taxonomy.qza

qiime metadata tabulate \
  --m-input-file step.04.taxonomy.qza \
  --o-visualization step.04.taxonomy.qzv

qiime taxa barplot \
  --i-table step.02.a.deblur.table.qza \
  --i-taxonomy step.04.taxonomy.qza \
  --m-metadata-file "$METADATA" \
  --o-visualization step.04.taxonomy.barplot.qzv

## STEP 5 — Export feature table to BIOM/TSV

module load qiime2/2024.10

qiime tools export \
  --input-path step.02.a.deblur.table.qza \
  --output-path biomresults

biom convert -i biomresults/feature-table.biom \
  -o biomresults/feature-table.txt \
  --to-tsv

done

# STEP 5 — Export feature table to BIOM/TSV
module load qiime2/2024.10

qiime tools export \
  --input-path step.02.a.deblur.table.qza \
  --output-path biomresults

biom convert -i biomresults/feature-table.biom \
  -o biomresults/feature-table.txt \
  --to-tsv

date

## STEP 6 — (Optional) Metadata pointers / version note

# This script loops through a set of metadata columns and automatically filters samples 
# from a QIIME 2 feature table based on each unique value in those columns.
# For each metadata value, it:
#   1. Extracts unique entries from the metadata file.
#   2. Cleans and sanitizes the value to use in file naming.
#   3. Runs 'qiime feature-table filter-samples' to generate a filtered table 
#      containing only samples matching that metadata value.
# The filtered tables are saved in the ANCOM directory with descriptive filenames.

module load qiime2/2024.10

metadata="path/to/metadata"
columns_to_check=("sample_type" "diet" "sample_id" "habitat" "country")


for column in "${columns_to_check[@]}"
do
    # Get unique values in the column (skip header)
    values=$(tail -n +2 "$metadata" | cut -f"$(head -1 $metadata | tr '\t' '\n' | grep -n -w "$column" | cut -d: -f1)" | sort | uniq)
    
    for value in $values
    do
        # Escape single quotes in value
        safe_value=$(echo "$value" | sed "s/'/''/g")

        # Create a sanitized output name
        clean_value=$(echo "$value" | tr ' ' '_' | tr -cd '[:alnum:]_')

        qiime feature-table filter-samples \
            --i-table step.02.a.deblur.table.qza \
            --m-metadata-file $metadata \
            --p-where "[${column}] = '${safe_value}'" \
            --o-filtered-table "/blue/cmavian/thongthum.t/bat16s/ANCOM/${column}_${clean_value}_ancom.qza"

    done
done

date

