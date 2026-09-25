#!/usr/bin/env bash
# ==============================================================================
# Bacterial 16S rRNA (V3–V4) amplicon analysis
# — Cutadapt + QIIME 2 + DADA2 + GSR-DB taxonomy (VSEARCH consensus)
# ==============================================================================
# Script  : 04_16s_amplicon_pipeline.sh
# Authors : <AUTHOR NAMES>
# Paper   : <CITATION / DOI OF PUBLISHED ARTICLE>
# Licence : <LICENCE, e.g. MIT>
#
# Steps   : primer trimming (Cutadapt) -> import -> DADA2 ASVs -> phylogenetic
#           tree -> GSR-DB taxonomy (VSEARCH) -> alpha/beta diversity
#           -> figures and tables (04_16s_figures.py)
#
# Requirements (tested with QIIME 2 amplicon 2024.10.1 on Ubuntu / WSL2):
#   - QIIME 2 amplicon distribution, installed with:
#       conda env create -n qiime2-amplicon-2024.10 \
#         --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
#       conda activate qiime2-amplicon-2024.10
#   - Cutadapt inside that environment:
#       conda install -n qiime2-amplicon-2024.10 -c conda-forge -c bioconda cutadapt
#   - GSR-DB V3–V4 DADA2-formatted reference FASTA
#     (GSR-DB_V3-V4_cluster-1_seqs_dada2.fasta)
#
# Input:
#   - A folder of demultiplexed paired-end FASTQ files (one R1 + one R2 per sample)
#   - A tab-separated sample metadata file. The first column must be named
#     "sample-id" and its values must match the FASTQ file names with the
#     R1/R2 suffix removed. Do not name any other column "SampleID".
#
# Usage:
#   Edit Section 0, activate the QIIME 2 environment, then run:
#     bash 04_16s_amplicon_pipeline.sh
#   Steps whose output already exists are skipped, so an interrupted run can be
#   resumed by running the script again.
#
# Methods and references:
#   Primer trimming  Cutadapt, 5' 341F/805R, min length 100 nt
#                    (Martin 2011; Klindworth et al. 2013)
#   ASVs             DADA2 denoise-paired, trunc-len 0, trunc-q 2, maxEE 2
#                    (Callahan et al. 2016)
#   Workflow         QIIME 2 amplicon 2024.10 (Bolyen et al. 2019)
#   Taxonomy         GSR-DB V3–V4 cluster-1; VSEARCH 80% identity, 51% consensus
#                    (Molano et al. 2024, mSystems, doi:10.1128/msystems.00950-23;
#                    Yarza et al. 2014)
#   Phylogeny        MAFFT + FastTree (q2-phylogeny)
#   Alpha diversity  Shannon, Gini–Simpson, Faith PD, observed ASVs (Faith 1992)
#   Beta diversity   Weighted UniFrac and Bray–Curtis PCoA; Bray–Curtis UPGMA
#                    (Lozupone & Knight 2005)
#   Rarefaction      Minimum sample depth after DADA2 (all samples retained)
# ==============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# ── 0. CONFIGURATION — replace the placeholders before running ───────────────

FASTQ_DIR="<PATH_TO_RAW_FASTQ_FOLDER>"
METADATA="<PATH_TO_SAMPLE_METADATA>.tsv"
WORK_DIR="<PATH_TO_WORKING_FOLDER>"
GSR_FASTA="<PATH_TO_GSR_DB_V3-V4_DADA2_FASTA>"

# File-name endings of the forward and reverse reads, e.g. "_R1.fastq.gz"
R1_SUFFIX="<R1_FILE_SUFFIX>"
R2_SUFFIX="<R2_FILE_SUFFIX>"

# Metadata columns used by the figure script (see 04_16s_figures.py)
LABEL_COLUMN="<METADATA_LABEL_COLUMN>"     # sample labels on figure axes
GROUP_COLUMN="<METADATA_GROUP_COLUMN>"     # colours/markers in figures
# Optional: genera always shown first in the heatmap, comma-separated
FOCUS_GENERA=""

THREADS=4

# Primers: 341F / 805R (Klindworth et al. 2013)
FWD_PRIMER="CCTACGGGNGGCWGCAG"
REV_PRIMER="GACTACHVGGGTATCTAATCC"
MIN_LENGTH=100

# DADA2 (primers removed first; quality-based truncation)
TRUNC_LEN_F=0
TRUNC_LEN_R=0
TRUNC_Q=2
MAX_EE_F=2
MAX_EE_R=2

# Taxonomy (VSEARCH consensus)
PERC_IDENTITY=0.80
MIN_CONSENSUS=0.51
MAXACCEPTS=10

# Rarefaction depth: leave empty to use the minimum sample depth
# (all samples retained), or set an integer.
SAMPLING_DEPTH=""


# ── Checks and setup ─────────────────────────────────────────────────────────

for v in FASTQ_DIR METADATA WORK_DIR GSR_FASTA R1_SUFFIX R2_SUFFIX; do
  if [[ "${!v}" == *"<"* ]]; then
    echo "ERROR: set $v in Section 0 before running." >&2; exit 1
  fi
done
command -v qiime    >/dev/null || { echo "ERROR: qiime not found — activate the QIIME 2 environment." >&2; exit 1; }
command -v cutadapt >/dev/null || { echo "ERROR: cutadapt not found in the QIIME 2 environment." >&2; exit 1; }
[[ -d "$FASTQ_DIR" ]] || { echo "ERROR: FASTQ_DIR not found: $FASTQ_DIR" >&2; exit 1; }
[[ -f "$METADATA"  ]] || { echo "ERROR: METADATA not found: $METADATA" >&2; exit 1; }
[[ -f "$GSR_FASTA" ]] || { echo "ERROR: GSR_FASTA not found: $GSR_FASTA" >&2; exit 1; }

header=$(head -n 1 "$METADATA")
[[ "$(cut -f1 <<< "$header")" == "sample-id" ]] || { echo "ERROR: first metadata column must be 'sample-id'." >&2; exit 1; }
if tr '\t' '\n' <<< "$header" | grep -qx "SampleID"; then
  echo "ERROR: 'SampleID' is a reserved name in QIIME 2 — rename that metadata column." >&2; exit 1
fi

METADATA=$(realpath "$METADATA")
FASTQ_DIR=$(realpath "$FASTQ_DIR")
GSR_FASTA=$(realpath "$GSR_FASTA")
mkdir -p "$WORK_DIR"/{qiime2,trimmed,results/{png,csv,qzv}}
cd "$WORK_DIR"
Q=qiime2
CM=$Q/core-metrics

export TMPDIR="$WORK_DIR/tmp"; export TEMP="$TMPDIR"; export TMP="$TMPDIR"
mkdir -p "$TMPDIR"

# Run a step only if its output does not exist yet
step() {
  local out="$1" name="$2"; shift 2
  if [[ -e "$out" ]]; then echo "-- skip  $name (exists: $out)"; return 0; fi
  echo "-- run   $name"; "$@"
}

echo "=== 16S pipeline started: $(date) ==="
qiime --version | head -n 1
echo "cutadapt $(cutadapt --version)"


# ── 1. Primer trimming (Cutadapt) and manifest of trimmed reads ──────────────

trim_reads() {
  python3 - "$FASTQ_DIR" "$R1_SUFFIX" "$R2_SUFFIX" trimmed "$Q/manifest_trimmed.tsv" \
            "$FWD_PRIMER" "$REV_PRIMER" "$MIN_LENGTH" "$THREADS" <<'PY'
import glob, os, subprocess, sys
fq_dir, r1, r2, trim, manifest, fwd, rev, min_len, threads = sys.argv[1:]
files = sorted(glob.glob(os.path.join(fq_dir, "*" + r1)))
if not files:
    sys.exit(f"ERROR: no files ending in {r1} in {fq_dir}")
rows = []
for f in files:
    sid = os.path.basename(f)[: -len(r1)]
    r_in = f[: -len(r1)] + r2
    if not os.path.exists(r_in):
        sys.exit(f"ERROR: reverse read missing for {sid}: {r_in}")
    f_out = os.path.join(trim, sid + "_R1.fastq.gz")
    r_out = os.path.join(trim, sid + "_R2.fastq.gz")
    if not (os.path.exists(f_out) and os.path.exists(r_out) and os.path.getsize(f_out) > 1000):
        print("cutadapt", sid, flush=True)
        subprocess.check_call(["cutadapt", "-j", threads, "-g", fwd, "-G", rev,
                               "--minimum-length", min_len,
                               "-o", f_out, "-p", r_out, f, r_in])
    rows.append((sid, os.path.abspath(f_out), os.path.abspath(r_out)))
with open(manifest, "w") as fh:
    fh.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
    for sid, a, b in rows:
        fh.write(f"{sid}\t{a}\t{b}\n")
print("Trimmed samples:", len(rows))
PY
}
step "$Q/manifest_trimmed.tsv" "Cutadapt" trim_reads


# ── 2. Import trimmed reads ──────────────────────────────────────────────────

step "$Q/demux-trimmed.qza" "import" \
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$Q/manifest_trimmed.tsv" \
    --output-path "$Q/demux-trimmed.qza" \
    --input-format PairedEndFastqManifestPhred33V2

step "results/qzv/demux-trimmed.qzv" "demux summary" \
  qiime demux summarize --i-data "$Q/demux-trimmed.qza" \
    --o-visualization results/qzv/demux-trimmed.qzv


# ── 3. DADA2 denoising -> ASVs ───────────────────────────────────────────────

step "$Q/table.qza" "DADA2" \
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$Q/demux-trimmed.qza" \
    --p-trunc-len-f "$TRUNC_LEN_F" \
    --p-trunc-len-r "$TRUNC_LEN_R" \
    --p-trunc-q "$TRUNC_Q" \
    --p-max-ee-f "$MAX_EE_F" \
    --p-max-ee-r "$MAX_EE_R" \
    --p-n-threads "$THREADS" \
    --o-table "$Q/table.qza" \
    --o-representative-sequences "$Q/rep-seqs.qza" \
    --o-denoising-stats "$Q/dada2-stats.qza"

step "results/qzv/dada2-stats.qzv" "DADA2 stats" \
  qiime metadata tabulate --m-input-file "$Q/dada2-stats.qza" \
    --o-visualization results/qzv/dada2-stats.qzv
step "results/qzv/table.qzv" "table summary" \
  qiime feature-table summarize --i-table "$Q/table.qza" \
    --m-sample-metadata-file "$METADATA" --o-visualization results/qzv/table.qzv


# ── 4. Phylogenetic tree (MAFFT + FastTree) ──────────────────────────────────

step "$Q/rooted-tree.qza" "phylogeny" \
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$Q/rep-seqs.qza" \
    --p-n-threads "$THREADS" \
    --o-alignment "$Q/aligned-rep-seqs.qza" \
    --o-masked-alignment "$Q/masked-aligned-rep-seqs.qza" \
    --o-tree "$Q/unrooted-tree.qza" \
    --o-rooted-tree "$Q/rooted-tree.qza"


# ── 5. GSR-DB reference (unique IDs) and taxonomy ────────────────────────────
# GSR-DB DADA2 FASTA headers carry the taxonomy and are not unique, so each
# reference is renamed gsr_000001, gsr_000002, ... with a separate taxonomy table.

format_gsr() {
  python3 - "$GSR_FASTA" "$Q/gsr-ref-seqs.fasta" "$Q/gsr-ref-taxonomy.tsv" <<'PY'
import sys
src, seq_out, tax_out = sys.argv[1:]
PREFIXES = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]

def parse_tax(h):
    h = h.strip().lstrip(">")
    if "tax=" in h:
        h = h.split("tax=", 1)[1]
    elif "\t" in h:
        h = h.split("\t", 1)[1]
    h = h.replace(",", ";")
    parts = [p.strip() for p in h.split(";") if p.strip()]
    out = [p if "__" in p else (PREFIXES[i] + p if i < len(PREFIXES) else p)
           for i, p in enumerate(parts)]
    return ";".join(out) if out else "Unassigned"

n, header, chunks = 0, None, []
with open(src) as fh, open(seq_out, "w") as so, open(tax_out, "w") as to:
    to.write("Feature ID\tTaxon\n")
    def flush():
        global n
        if header is None:
            return
        n += 1
        seq = "".join(chunks).replace(" ", "").replace("U", "T")
        if seq:
            sid = "gsr_%06d" % n
            so.write(f">{sid}\n{seq}\n")
            to.write(f"{sid}\t{parse_tax(header)}\n")
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            header, chunks = line[1:], []
        else:
            chunks.append(line.strip())
    flush()
print("GSR-DB references:", n)
PY
}
step "$Q/gsr-ref-taxonomy.tsv" "format GSR-DB" format_gsr

step "$Q/gsr-ref-seqs.qza" "import GSR-DB sequences" \
  qiime tools import --type 'FeatureData[Sequence]' \
    --input-path "$Q/gsr-ref-seqs.fasta" --output-path "$Q/gsr-ref-seqs.qza"
step "$Q/gsr-ref-taxonomy.qza" "import GSR-DB taxonomy" \
  qiime tools import --type 'FeatureData[Taxonomy]' --input-format TSVTaxonomyFormat \
    --input-path "$Q/gsr-ref-taxonomy.tsv" --output-path "$Q/gsr-ref-taxonomy.qza"

step "$Q/taxonomy.qza" "classification" \
  qiime feature-classifier classify-consensus-vsearch \
    --i-query "$Q/rep-seqs.qza" \
    --i-reference-reads "$Q/gsr-ref-seqs.qza" \
    --i-reference-taxonomy "$Q/gsr-ref-taxonomy.qza" \
    --p-maxaccepts "$MAXACCEPTS" \
    --p-perc-identity "$PERC_IDENTITY" \
    --p-min-consensus "$MIN_CONSENSUS" \
    --p-threads "$THREADS" \
    --o-classification "$Q/taxonomy.qza" \
    --o-search-results "$Q/taxonomy-hits.qza"


# ── 6. Diversity ─────────────────────────────────────────────────────────────

if [[ -z "$SAMPLING_DEPTH" ]]; then
  # Only the integer goes to stdout; the summary line goes to stderr
  SAMPLING_DEPTH=$(python3 - "$Q/table.qza" <<'PY'
import sys
import pandas as pd
from qiime2 import Artifact
depths = Artifact.load(sys.argv[1]).view(pd.DataFrame).sum(axis=1).astype(float)
print(f"   sample depths min/median/max: {depths.min():.0f} / {depths.median():.0f} / {depths.max():.0f}",
      file=sys.stderr)
print(int(depths.min()))
PY
)
fi
[[ "$SAMPLING_DEPTH" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid sampling depth: $SAMPLING_DEPTH" >&2; exit 1; }
echo "   rarefaction depth: $SAMPLING_DEPTH"

step "$CM" "core metrics (depth $SAMPLING_DEPTH)" \
  qiime diversity core-metrics-phylogenetic \
    --i-table "$Q/table.qza" \
    --i-phylogeny "$Q/rooted-tree.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA" \
    --output-dir "$CM"

step "$Q/shannon.qza" "Shannon" \
  qiime diversity alpha --i-table "$CM/rarefied_table.qza" \
    --p-metric shannon --o-alpha-diversity "$Q/shannon.qza"
step "$Q/simpson.qza" "Simpson" \
  qiime diversity alpha --i-table "$CM/rarefied_table.qza" \
    --p-metric simpson --o-alpha-diversity "$Q/simpson.qza"
step "$Q/faith_pd.qza" "Faith PD" \
  qiime diversity alpha-phylogenetic --i-table "$CM/rarefied_table.qza" \
    --i-phylogeny "$Q/rooted-tree.qza" \
    --p-metric faith_pd --o-alpha-diversity "$Q/faith_pd.qza"


# ── 7. Figures and tables ────────────────────────────────────────────────────

python3 "$SCRIPT_DIR/04_16s_figures.py" \
  --qiime-dir "$WORK_DIR/$Q" --out-dir "$WORK_DIR/results" \
  --metadata "$METADATA" --label-column "$LABEL_COLUMN" --group-column "$GROUP_COLUMN" \
  --focus-genera "$FOCUS_GENERA"

echo "=== 16S pipeline complete: $(date) ==="
echo "    QIIME 2 artefacts: $WORK_DIR/$Q   |   figures/tables: $WORK_DIR/results"
