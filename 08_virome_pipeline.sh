#!/usr/bin/env bash
#==============================================================================
# 08_virome_pipeline.sh
# Plant virome from sequencing reads: quality/adapter trimming (fastp) ->
# taxonomic classification (Kraken2, RefSeq viral) -> species count matrix ->
# diversity tables and figures (08_virome_analysis.R)
#
# Authors : <AUTHOR NAMES>
# Paper   : <CITATION / DOI OF PUBLISHED ARTICLE>
# Licence : <LICENCE, e.g. MIT>
#
# Requirements: run 08_virome_setup.sh once first.
#
# Input:
#   - A folder of FASTQ files, one file per sample (single-end or pre-merged
#     reads), named <sample-id>.fastq.gz
#   - Optional sample metadata (tab-separated): first column "sample-id"
#     (file name without .fastq.gz), plus optional columns "Label" (name shown
#     in tables and figures) and "Group" (e.g. symptomatic/asymptomatic). Row
#     order sets the sample order in tables and figures.
#
# Usage:  bash 08_virome_pipeline.sh
#   Steps whose output already exists are skipped, so an interrupted run can be
#   resumed by running the script again.
#
# Software to cite: fastp (Chen et al. 2018); Kraken2 (Wood et al. 2019);
#   vegan (Oksanen et al.); ape (Paradis & Schliep 2019).
#==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# -------------------- SETTINGS (edit) --------------------
INPUT_DIR="<PATH_TO_FASTQ_FOLDER>"
OUTDIR="<PATH_TO_OUTPUT_FOLDER>"
KRAKEN_DB="<PATH_TO_KRAKEN2_VIRAL_DB_FOLDER>"
METADATA=""                                   # optional: path to sample metadata TSV
THREADS=4
ENV_NAME="potato_virome"

# fastp
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"   # Illumina TruSeq read 1 adapter
MIN_QUAL=20
MIN_LENGTH=50

# -------------------- ENVIRONMENT --------------------
for v in INPUT_DIR OUTDIR KRAKEN_DB; do
  if [[ "${!v}" == *"<"* ]]; then echo "ERROR: set $v at the top of this script." >&2; exit 1; fi
done

if ! command -v kraken2 >/dev/null 2>&1; then
  # shellcheck disable=SC1090,SC1091
  for c in "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh"; do
    [ -f "$c" ] && source "$c" && break
  done
  if conda env list 2>/dev/null | grep -qE "^${ENV_NAME}\s"; then
    conda activate "${ENV_NAME}"
  else
    echo "ERROR: conda environment '${ENV_NAME}' not found. Run 08_virome_setup.sh first." >&2; exit 1
  fi
fi

# -------------------- CHECKS --------------------
[ -d "${INPUT_DIR}" ] || { echo "ERROR: input folder not found: ${INPUT_DIR}" >&2; exit 1; }
[ -f "${KRAKEN_DB}/hash.k2d" ] || { echo "ERROR: Kraken2 database files missing in: ${KRAKEN_DB}" >&2; exit 1; }
if [ -n "${METADATA}" ]; then
  [ -f "${METADATA}" ] || { echo "ERROR: METADATA not found: ${METADATA}" >&2; exit 1; }
  [[ "$(head -n 1 "${METADATA}" | cut -f1)" == "sample-id" ]] || \
    { echo "ERROR: first metadata column must be 'sample-id'." >&2; exit 1; }
fi

mkdir -p "${OUTDIR}/01_trimmed" "${OUTDIR}/02_kraken" "${OUTDIR}/00_logs"
LOG="${OUTDIR}/00_logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

echo "============================================================"
echo " Virome pipeline"
echo " Started : $(date)"
echo " Input   : ${INPUT_DIR}"
echo " Output  : ${OUTDIR}"
echo " Database: ${KRAKEN_DB} ($(cat "${KRAKEN_DB}/DATABASE_VERSION.txt" 2>/dev/null || echo "version not recorded"))"
echo "============================================================"

# -------------------- STEP 1: FASTP TRIMMING --------------------
echo ""
echo ">>> [1/4] Adapter and quality trimming (fastp)"
shopt -s nullglob
raw_fastqs=("${INPUT_DIR}"/*.fastq.gz)
shopt -u nullglob
[ ${#raw_fastqs[@]} -gt 0 ] || { echo "ERROR: no .fastq.gz files in ${INPUT_DIR}" >&2; exit 1; }

for fq in "${raw_fastqs[@]}"; do
  sample=$(basename "${fq}" .fastq.gz)
  out_fq="${OUTDIR}/01_trimmed/${sample}.trimmed.fastq.gz"
  if [ -s "${out_fq}" ]; then echo "  [skip] ${sample}"; continue; fi
  echo "  [run]  fastp ${sample}"
  fastp \
    -i "${fq}" \
    -o "${out_fq}" \
    -h "${OUTDIR}/01_trimmed/${sample}.fastp.html" \
    -j "${OUTDIR}/01_trimmed/${sample}.fastp.json" \
    --adapter_sequence "${ADAPTER}" \
    --qualified_quality_phred "${MIN_QUAL}" \
    --length_required "${MIN_LENGTH}" \
    --trim_poly_g \
    --thread "${THREADS}" \
    --compression 6
done

# -------------------- STEP 2: KRAKEN2 CLASSIFICATION --------------------
echo ""
echo ">>> [2/4] Taxonomic classification (Kraken2)"
for trim in "${OUTDIR}/01_trimmed/"*.trimmed.fastq.gz; do
  sample=$(basename "${trim}" .trimmed.fastq.gz)
  report="${OUTDIR}/02_kraken/${sample}.kraken2.report"
  if [ -s "${report}" ]; then echo "  [skip] ${sample}"; continue; fi
  echo "  [run]  kraken2 ${sample}"
  kraken2 \
    --db "${KRAKEN_DB}" \
    --threads "${THREADS}" \
    --report "${report}" \
    --output "${OUTDIR}/02_kraken/${sample}.kraken2.out" \
    "${trim}"
done

# -------------------- STEP 3: SPECIES COUNT MATRIX --------------------
# One row per species (Kraken2 rank code "S"), using the clade read count,
# which already includes reads assigned to strains/sub-species (S1, S2, ...).
# Sub-species rows are not added again, so no read is counted twice.
echo ""
echo ">>> [3/4] Building species count matrix"
export OUTDIR
python3 - <<'PY'
import glob, os
from collections import defaultdict

kraken_dir = os.path.join(os.environ["OUTDIR"], "02_kraken")
reports = sorted(glob.glob(os.path.join(kraken_dir, "*.kraken2.report")))
if not reports:
    raise SystemExit(f"ERROR: no Kraken2 reports in {kraken_dir}")

mat, samples = defaultdict(dict), []
for rep in reports:
    sample = os.path.basename(rep).replace(".kraken2.report", "")
    samples.append(sample)
    with open(rep) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or parts[3].strip() != "S":
                continue
            try:
                clade_reads = int(parts[1].strip())
            except ValueError:
                continue
            if clade_reads > 0:
                name = parts[5].strip()
                mat[name][sample] = mat[name].get(sample, 0) + clade_reads

samples = sorted(samples)
out_tsv = os.path.join(kraken_dir, "species_abundance_counts_raw.tsv")
with open(out_tsv, "w") as out:
    out.write("taxon\t" + "\t".join(samples) + "\n")
    for t in sorted(mat):
        out.write("\t".join([t] + [str(mat[t].get(s, 0)) for s in samples]) + "\n")
print(f"Matrix: {out_tsv} ({len(mat)} species x {len(samples)} samples)")
PY

# -------------------- STEP 4: DIVERSITY, TABLES AND FIGURES --------------------
echo ""
echo ">>> [4/4] Diversity analysis, tables and figures (R)"
export METADATA
Rscript "${SCRIPT_DIR}/08_virome_analysis.R"

echo ""
echo "============================================================"
echo " Pipeline complete: $(date)"
echo " Log: ${LOG}"
echo "============================================================"
