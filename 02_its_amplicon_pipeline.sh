#!/usr/bin/env bash
# ==============================================================================
# Fungal ITS amplicon analysis — QIIME 2 + DADA2 + UNITE (VSEARCH consensus)
# ==============================================================================
# Script  : 02_its_amplicon_pipeline.sh
# Authors : <AUTHOR NAMES>
# Paper   : <CITATION / DOI OF PUBLISHED ARTICLE>
# Licence : <LICENCE, e.g. MIT>
#
# Steps   : import paired-end reads -> (optional) primer trimming -> DADA2 ASVs
#           -> phylogenetic tree -> UNITE taxonomy (VSEARCH) -> alpha/beta
#           diversity -> figures/tables -> ASV sequence export
#
# Requirements (tested with QIIME 2 amplicon 2024.10.1 on Ubuntu / WSL2):
#   - 4+ CPU cores, >= 16 GB RAM, ~40 GB free disk
#   - QIIME 2 amplicon distribution, installed with:
#       conda env create -n qiime2-amplicon-2024.10 \
#         --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
#       conda activate qiime2-amplicon-2024.10
#   - BLAST+ only if RUN_BLAST=true:  conda install -c bioconda blast
#
# Input:
#   - A folder of demultiplexed paired-end FASTQ files (one R1 + one R2 per sample)
#   - A tab-separated sample metadata file. The first column must be named
#     "sample-id" and its values must match the FASTQ file names with the
#     R1/R2 suffix removed. Do not name any other column "SampleID".
#     See sample_metadata_template.tsv.
#
# Usage:
#   Edit Section 0, activate the QIIME 2 environment, then run:
#     bash 02_its_amplicon_pipeline.sh
#   Steps whose output already exists are skipped, so an interrupted run can be
#   resumed by running the script again.
#
# Software to cite: QIIME 2 (Bolyen et al. 2019); DADA2 (Callahan et al. 2016);
#   Cutadapt (Martin 2011) if used; UNITE (Abarenkov et al.); VSEARCH
#   (Rognes et al. 2016) via q2-feature-classifier; MAFFT; FastTree.
# ==============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# ── 0. CONFIGURATION — replace the placeholders before running ───────────────

FASTQ_DIR="<PATH_TO_RAW_FASTQ_FOLDER>"
METADATA="<PATH_TO_SAMPLE_METADATA>.tsv"
WORK_DIR="<PATH_TO_WORKING_FOLDER>"

# File-name endings of the forward and reverse reads, e.g. "_R1.fastq.gz"
R1_SUFFIX="<R1_FILE_SUFFIX>"
R2_SUFFIX="<R2_FILE_SUFFIX>"

# Metadata column used for group comparisons (beta-group-significance)
GROUP_COLUMN="<METADATA_GROUP_COLUMN>"

# Compute resources
THREADS=4
CLASSIFY_THREADS=2

# Primer trimming (Cutadapt). Primers: ITS1F (forward) / ITS2 (reverse)
RUN_CUTADAPT=false
PRIMER_F="CTTGGTCATTTAGAGGAAGTAA"
PRIMER_R="GCTGCGTTCTTCATCGATGC"

# DADA2 — choose truncation lengths from the quality plots in demux.qzv
TRUNC_LEN_F=250
TRUNC_LEN_R=200
TRIM_LEFT_F=0
TRIM_LEFT_R=0
MAX_EE_F=2
MAX_EE_R=2

# Taxonomy — UNITE QIIME release v10.0 (04.04.2024), dynamic, developer files
UNITE_URL="https://s3.hpc.ut.ee/plutof-public/original/db1d6ddb-a35d-48c5-8b1a-ad9dd3310c6d.tgz"
UNITE_FASTA="developer/sh_refs_qiime_ver10_dynamic_04.04.2024_dev.fasta"
UNITE_TAX="developer/sh_taxonomy_qiime_ver10_dynamic_04.04.2024_dev.txt"
PERC_IDENTITY=0.97

# Rarefaction depth — choose from the table summary so that samples are retained
SAMPLING_DEPTH=3000

# Optional: export and BLAST the ASVs of one genus (leave empty to skip)
EXPORT_GENUS=""
RUN_BLAST=false


# ── Checks and setup ─────────────────────────────────────────────────────────

for v in FASTQ_DIR METADATA WORK_DIR R1_SUFFIX R2_SUFFIX GROUP_COLUMN; do
  if [[ "${!v}" == *"<"* ]]; then
    echo "ERROR: set $v in Section 0 before running." >&2; exit 1
  fi
done
command -v qiime >/dev/null || { echo "ERROR: qiime not found — activate the QIIME 2 environment." >&2; exit 1; }
[[ -d "$FASTQ_DIR" ]] || { echo "ERROR: FASTQ_DIR not found: $FASTQ_DIR" >&2; exit 1; }
[[ -f "$METADATA"  ]] || { echo "ERROR: METADATA not found: $METADATA" >&2; exit 1; }

header=$(head -n 1 "$METADATA")
[[ "$(cut -f1 <<< "$header")" == "sample-id" ]] || { echo "ERROR: first metadata column must be 'sample-id'." >&2; exit 1; }
tr '\t' '\n' <<< "$header" | grep -qx "$GROUP_COLUMN" || { echo "ERROR: column '$GROUP_COLUMN' not in metadata." >&2; exit 1; }
if tr '\t' '\n' <<< "$header" | grep -qx "SampleID"; then
  echo "ERROR: 'SampleID' is a reserved name in QIIME 2 — rename that metadata column." >&2; exit 1
fi

METADATA=$(realpath "$METADATA")
FASTQ_DIR=$(realpath "$FASTQ_DIR")
mkdir -p "$WORK_DIR"/{qiime2,db,results/{csv,figures,all_fungi}}
cd "$WORK_DIR"
Q=qiime2
CM=$Q/core-metrics-results

# A dedicated temp folder avoids I/O errors from a small /tmp (e.g. under WSL)
export TMPDIR="$WORK_DIR/tmp"; export TEMP="$TMPDIR"; export TMP="$TMPDIR"
mkdir -p "$TMPDIR"

# Run a step only if its output does not exist yet
step() {
  local out="$1" name="$2"; shift 2
  if [[ -e "$out" ]]; then echo "-- skip  $name (exists: $out)"; return 0; fi
  echo "-- run   $name"; "$@"
}

echo "=== ITS pipeline started: $(date) ==="


# ── 1. Import reads ──────────────────────────────────────────────────────────

make_manifest() {
  python3 - "$FASTQ_DIR" "$R1_SUFFIX" "$R2_SUFFIX" "$Q/manifest.tsv" <<'PY'
import glob, os, sys
fq_dir, r1, r2, out = sys.argv[1:]
files = sorted(glob.glob(os.path.join(fq_dir, "*" + r1)))
if not files:
    sys.exit(f"ERROR: no files ending in {r1} in {fq_dir}")
with open(out, "w") as fh:
    fh.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
    for f in files:
        sid = os.path.basename(f)[: -len(r1)]
        rev = f[: -len(r1)] + r2
        if not os.path.exists(rev):
            sys.exit(f"ERROR: reverse read missing for {sid}: {rev}")
        fh.write(f"{sid}\t{os.path.abspath(f)}\t{os.path.abspath(rev)}\n")
print(f"Manifest created with {len(files)} samples")
PY
}
step "$Q/manifest.tsv" "manifest" make_manifest

step "$Q/demux.qza" "import" \
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$Q/manifest.tsv" \
    --output-path "$Q/demux.qza" \
    --input-format PairedEndFastqManifestPhred33V2

step "$Q/demux.qzv" "demux summary" \
  qiime demux summarize --i-data "$Q/demux.qza" --o-visualization "$Q/demux.qzv"


# ── 2. Primer trimming (optional) ────────────────────────────────────────────

DADA2_INPUT="$Q/demux.qza"
if [[ "$RUN_CUTADAPT" == true ]]; then
  if step "$Q/demux-trimmed.qza" "cutadapt" \
       qiime cutadapt trim-paired \
         --i-demultiplexed-sequences "$Q/demux.qza" \
         --p-front-f "$PRIMER_F" \
         --p-front-r "$PRIMER_R" \
         --p-error-rate 0.1 \
         --p-cores "$THREADS" \
         --o-trimmed-sequences "$Q/demux-trimmed.qza"; then
    DADA2_INPUT="$Q/demux-trimmed.qza"
  else
    echo "WARNING: Cutadapt failed — DADA2 will run on untrimmed reads." >&2
  fi
fi
echo "   DADA2 input: $DADA2_INPUT"


# ── 3. DADA2 denoising -> ASVs ───────────────────────────────────────────────

step "$Q/table.qza" "DADA2" \
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$DADA2_INPUT" \
    --p-trunc-len-f "$TRUNC_LEN_F" \
    --p-trunc-len-r "$TRUNC_LEN_R" \
    --p-trim-left-f "$TRIM_LEFT_F" \
    --p-trim-left-r "$TRIM_LEFT_R" \
    --p-max-ee-f "$MAX_EE_F" \
    --p-max-ee-r "$MAX_EE_R" \
    --p-n-threads "$THREADS" \
    --o-table "$Q/table.qza" \
    --o-representative-sequences "$Q/rep-seqs.qza" \
    --o-denoising-stats "$Q/denoising-stats.qza"

step "$Q/denoising-stats.qzv" "DADA2 stats" \
  qiime metadata tabulate --m-input-file "$Q/denoising-stats.qza" --o-visualization "$Q/denoising-stats.qzv"
step "$Q/table.qzv" "table summary" \
  qiime feature-table summarize --i-table "$Q/table.qza" --o-visualization "$Q/table.qzv" \
    --m-sample-metadata-file "$METADATA"
step "$Q/rep-seqs.qzv" "ASV sequences" \
  qiime feature-table tabulate-seqs --i-data "$Q/rep-seqs.qza" --o-visualization "$Q/rep-seqs.qzv"


# ── 4. Phylogenetic tree (for UniFrac and Faith's PD) ────────────────────────

step "$Q/rooted-tree.qza" "phylogeny" \
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$Q/rep-seqs.qza" \
    --o-alignment "$Q/aligned-rep-seqs.qza" \
    --o-masked-alignment "$Q/masked-aligned-rep-seqs.qza" \
    --o-tree "$Q/unrooted-tree.qza" \
    --o-rooted-tree "$Q/rooted-tree.qza" \
    --p-n-threads "$THREADS"


# ── 5. UNITE taxonomy (VSEARCH consensus) ────────────────────────────────────

get_unite() {
  mkdir -p db/unite_raw
  wget -O db/unite_raw/unite_qiime.tgz "$UNITE_URL"
  tar -xzf db/unite_raw/unite_qiime.tgz -C db/unite_raw
}
step "db/unite_raw/$UNITE_FASTA" "UNITE download" get_unite

step "db/unite-ref-seqs.qza" "UNITE import (sequences)" \
  qiime tools import --type 'FeatureData[Sequence]' \
    --input-path "db/unite_raw/$UNITE_FASTA" \
    --output-path db/unite-ref-seqs.qza
step "db/unite-ref-tax.qza" "UNITE import (taxonomy)" \
  qiime tools import --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path "db/unite_raw/$UNITE_TAX" \
    --output-path db/unite-ref-tax.qza

step "$Q/taxonomy.qza" "classification" \
  qiime feature-classifier classify-consensus-vsearch \
    --i-query "$Q/rep-seqs.qza" \
    --i-reference-reads db/unite-ref-seqs.qza \
    --i-reference-taxonomy db/unite-ref-tax.qza \
    --p-perc-identity "$PERC_IDENTITY" \
    --p-threads "$CLASSIFY_THREADS" \
    --o-classification "$Q/taxonomy.qza" \
    --o-search-results "$Q/taxonomy-search.qza"

step "$Q/taxonomy.qzv" "taxonomy table" \
  qiime metadata tabulate --m-input-file "$Q/taxonomy.qza" --o-visualization "$Q/taxonomy.qzv"
step "$Q/taxa-barplot.qzv" "taxa barplot" \
  qiime taxa barplot --i-table "$Q/table.qza" --i-taxonomy "$Q/taxonomy.qza" \
    --m-metadata-file "$METADATA" --o-visualization "$Q/taxa-barplot.qzv"


# ── 6. Diversity analyses ────────────────────────────────────────────────────

step "$Q/exported-table/feature-table.biom" "export table" \
  qiime tools export --input-path "$Q/table.qza" --output-path "$Q/exported-table"
if command -v biom >/dev/null; then
  biom summarize-table -i "$Q/exported-table/feature-table.biom" > results/csv/table_depth_summary.txt
  echo "   Read depth per sample: results/csv/table_depth_summary.txt (check SAMPLING_DEPTH=$SAMPLING_DEPTH)"
fi

step "$CM" "core metrics (depth $SAMPLING_DEPTH)" \
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "$Q/rooted-tree.qza" \
    --i-table "$Q/table.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA" \
    --output-dir "$CM" \
    --p-n-jobs-or-threads "$CLASSIFY_THREADS"

# Simpson is not produced by core-metrics
step "$CM/simpson_vector.qza" "Simpson" \
  qiime diversity alpha --i-table "$CM/rarefied_table.qza" \
    --p-metric simpson --o-alpha-diversity "$CM/simpson_vector.qza"

for metric in shannon simpson faith_pd; do
  step "$Q/${metric}-group-significance.qzv" "alpha significance ($metric)" \
    qiime diversity alpha-group-significance \
      --i-alpha-diversity "$CM/${metric}_vector.qza" \
      --m-metadata-file "$METADATA" \
      --o-visualization "$Q/${metric}-group-significance.qzv"
done

for dist in bray_curtis weighted_unifrac; do
  step "$Q/${dist}-${GROUP_COLUMN}-significance.qzv" "beta significance ($dist by $GROUP_COLUMN)" \
    qiime diversity beta-group-significance \
      --i-distance-matrix "$CM/${dist}_distance_matrix.qza" \
      --m-metadata-file "$METADATA" \
      --m-metadata-column "$GROUP_COLUMN" \
      --p-pairwise \
      --o-visualization "$Q/${dist}-${GROUP_COLUMN}-significance.qzv"
done


# ── 7. Figures and tables ────────────────────────────────────────────────────
# Genus heatmap, alpha diversity, weighted UniFrac PCoA, Bray–Curtis UPGMA.

FIG_SCRIPT="$SCRIPT_DIR/<FIGURE_SCRIPT_NAME>.py"
if [[ -f "$FIG_SCRIPT" ]]; then
  python3 "$FIG_SCRIPT"
else
  echo "-- skip  figures (figure script not found: $FIG_SCRIPT)"
fi


# ── 8. Export ASV sequences ──────────────────────────────────────────────────

EXPORT_SCRIPT="$SCRIPT_DIR/<EXPORT_SCRIPT_NAME>.sh"
if [[ -f "$EXPORT_SCRIPT" ]]; then
  bash "$EXPORT_SCRIPT"
else
  echo "-- skip  full ASV export (export script not found: $EXPORT_SCRIPT)"
fi

if [[ -n "$EXPORT_GENUS" ]]; then
  G="results/$EXPORT_GENUS"
  mkdir -p "$G"
  step "$G/${EXPORT_GENUS}-seqs.qza" "filter $EXPORT_GENUS ASVs" \
    qiime taxa filter-seqs --i-sequences "$Q/rep-seqs.qza" --i-taxonomy "$Q/taxonomy.qza" \
      --p-include "$EXPORT_GENUS" --o-filtered-sequences "$G/${EXPORT_GENUS}-seqs.qza"
  step "$G/seqs/dna-sequences.fasta" "export $EXPORT_GENUS FASTA" \
    qiime tools export --input-path "$G/${EXPORT_GENUS}-seqs.qza" --output-path "$G/seqs"

  # Optional remote NCBI BLAST of the selected genus (after UNITE classification)
  if [[ "$RUN_BLAST" == true ]]; then
    step "$G/${EXPORT_GENUS}_NCBI_blast.tsv" "BLAST $EXPORT_GENUS" \
      blastn -query "$G/seqs/dna-sequences.fasta" -db nt -remote -task blastn \
        -evalue 1e-20 -max_target_seqs 5 -max_hsps 1 \
        -outfmt "6 qseqid sacc staxids sscinames pident length qcovs evalue bitscore stitle" \
        -out "$G/${EXPORT_GENUS}_NCBI_blast.tsv"
  fi
fi

echo "=== ITS pipeline complete: $(date) ==="
echo "    QIIME 2 artefacts: $WORK_DIR/$Q   |   results: $WORK_DIR/results"
echo "    View .qzv files at https://view.qiime2.org"
