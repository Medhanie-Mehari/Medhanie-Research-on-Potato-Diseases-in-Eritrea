#!/usr/bin/env python3
# ==============================================================================
# Fungal ITS amplicons: DADA2 ASVs -> remote NCBI nt BLASTn of every ASV
# (resumable) -> best-hit annotation tables
# ==============================================================================
# Script  : 03_its_ncbi_nt_blast.py
# Authors : <AUTHOR NAMES>
# Paper   : <CITATION / DOI OF PUBLISHED ARTICLE>
# Licence : <LICENCE, e.g. MIT>
#
# Input:
#   - A folder of demultiplexed paired-end FASTQ files (one R1 + one R2 per sample)
#   - A tab-separated sample metadata file whose first column is "sample-id",
#     with values matching the FASTQ file names minus the R1/R2 suffix.
#
# Method: reads are imported into QIIME 2, optionally primer-trimmed with
# Cutadapt, and denoised with DADA2 into ASVs. Each ASV is BLASTed individually against the NCBI nt database with
# remote blastn, in order of decreasing total read count. Completed ASV IDs are
# recorded, so the job can be restarted after an interruption and continues
# where it stopped.
#
# Requirements: QIIME 2 amplicon 2024.10 environment (provides pandas and the
#   qiime2 Python API) and BLAST+ (conda install -c bioconda blast).
#   Run only one instance at a time — NCBI may block parallel remote jobs.
#
# Usage (inside the activated QIIME 2 environment):
#   python 03_its_ncbi_nt_blast.py          # all steps: ASVs, BLAST, tables
#   python 03_its_ncbi_nt_blast.py asv      # import + DADA2 only
#   python 03_its_ncbi_nt_blast.py blast    # BLAST only (resumable)
#   python 03_its_ncbi_nt_blast.py tables   # tables only (works on partial results)
#   Steps whose output already exists are skipped, so the script can be
#   re-run after an interruption.
# ==============================================================================

import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

# ── 0. CONFIGURATION — replace the placeholders before running ───────────────

FASTQ_DIR = Path("<PATH_TO_RAW_FASTQ_FOLDER>")
METADATA = Path("<PATH_TO_SAMPLE_METADATA>.tsv")
WORK_DIR = Path("<PATH_TO_WORKING_FOLDER>")
LABEL_COLUMN = "<METADATA_LABEL_COLUMN>"              # e.g. Location; used to label samples

# File-name endings of the forward and reverse reads, e.g. "_R1.fastq.gz"
R1_SUFFIX = "<R1_FILE_SUFFIX>"
R2_SUFFIX = "<R2_FILE_SUFFIX>"

THREADS = 4

# Primer trimming (Cutadapt). Primers: ITS1F (forward) / ITS2 (reverse)
RUN_CUTADAPT = False
PRIMER_F = "CTTGGTCATTTAGAGGAAGTAA"
PRIMER_R = "GCTGCGTTCTTCATCGATGC"

# DADA2 — choose truncation lengths from the quality plots in demux.qzv
TRUNC_LEN_F = 250
TRUNC_LEN_R = 200
TRIM_LEFT_F = 0
TRIM_LEFT_R = 0
MAX_EE_F = 2
MAX_EE_R = 2

Q_DIR = WORK_DIR / "qiime2"
MANIFEST = Q_DIR / "manifest.tsv"
DEMUX_QZA = Q_DIR / "demux.qza"
TRIMMED_QZA = Q_DIR / "demux-trimmed.qza"
TABLE_QZA = WORK_DIR / "qiime2" / "table.qza"
REPSEQS_QZA = WORK_DIR / "qiime2" / "rep-seqs.qza"
OUT_DIR = WORK_DIR / "results" / "ncbi_nt_blast"

# BLAST settings
BLAST_DB = "nt"
EVALUE = "1e-10"
MAX_TARGET_SEQS = 3
MAX_HSPS = 1
TIMEOUT_S = 1800          # per ASV (30 min)
ATTEMPTS = 2              # tries per ASV before it is recorded as failed
PAUSE_RETRY_S = 15
PAUSE_BETWEEN_S = 5

OUTFMT = "6 qseqid sacc staxids sscinames pident length qcovs evalue bitscore stitle"
HIT_COLS = ["FeatureID", "NCBI_acc", "TaxID", "SciName", "pident",
            "length", "qcovs", "evalue", "bitscore", "title"]

FASTA = OUT_DIR / "seqs" / "dna-sequences.fasta"
RAW = OUT_DIR / "nt_blast_raw.tsv"
DONE = OUT_DIR / "nt_blast_done_ids.txt"
FAILED = OUT_DIR / "nt_blast_failed_ids.txt"
ONE = OUT_DIR / "one.fasta"


# ── Helpers ──────────────────────────────────────────────────────────────────

def check_config():
    for name, val in [("FASTQ_DIR", FASTQ_DIR), ("METADATA", METADATA), ("WORK_DIR", WORK_DIR),
                      ("R1_SUFFIX", R1_SUFFIX), ("R2_SUFFIX", R2_SUFFIX)]:
        if "<" in str(val):
            sys.exit(f"ERROR: set {name} in Section 0 before running.")
    if not FASTQ_DIR.is_dir():
        sys.exit(f"ERROR: FASTQ_DIR not found: {FASTQ_DIR}")
    if not METADATA.is_file():
        sys.exit(f"ERROR: METADATA not found: {METADATA}")
    header = METADATA.read_text().splitlines()[0].split("\t")
    if header[0] != "sample-id":
        sys.exit("ERROR: first metadata column must be 'sample-id'.")
    Q_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)


def step(output, name, cmd):
    """Run a command only if its output does not exist yet."""
    if Path(output).exists():
        print(f"-- skip  {name} (exists: {output})", flush=True)
        return True
    print(f"-- run   {name}", flush=True)
    return subprocess.run([str(c) for c in cmd]).returncode == 0


def read_fasta(path):
    seqs, sid, bits = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if sid:
                seqs[sid] = "".join(bits)
            sid, bits = line[1:].split()[0], []
        else:
            bits.append(line.strip())
    if sid:
        seqs[sid] = "".join(bits)
    return seqs


def load_table(feature_ids):
    """ASV table as samples (rows) x features (columns), numeric."""
    from qiime2 import Artifact
    table = Artifact.load(str(TABLE_QZA)).view(pd.DataFrame)
    table.index = table.index.astype(str)
    table.columns = table.columns.astype(str)
    # Orient so that features are columns
    if len(set(table.index) & feature_ids) > len(set(table.columns) & feature_ids):
        table = table.T
    return table.apply(pd.to_numeric, errors="coerce").fillna(0)


def read_ids(path):
    return set(path.read_text().split()) if path.exists() else set()


def append_line(path, text):
    with path.open("a") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")


# ── 1. Import reads and denoise with DADA2 -> ASVs ───────────────────────────

def make_manifest():
    if MANIFEST.exists():
        print(f"-- skip  manifest (exists: {MANIFEST})")
        return
    r1_files = sorted(FASTQ_DIR.glob("*" + R1_SUFFIX))
    if not r1_files:
        sys.exit(f"ERROR: no files ending in {R1_SUFFIX} in {FASTQ_DIR}")
    lines = ["sample-id\tforward-absolute-filepath\treverse-absolute-filepath"]
    for f in r1_files:
        sid = f.name[: -len(R1_SUFFIX)]
        rev = f.with_name(sid + R2_SUFFIX)
        if not rev.exists():
            sys.exit(f"ERROR: reverse read missing for {sid}: {rev}")
        lines.append(f"{sid}\t{f.resolve()}\t{rev.resolve()}")
    MANIFEST.write_text("\n".join(lines) + "\n")
    print(f"Manifest created with {len(r1_files)} samples")


def make_asvs():
    make_manifest()
    ok = step(DEMUX_QZA, "import", [
        "qiime", "tools", "import",
        "--type", "SampleData[PairedEndSequencesWithQuality]",
        "--input-path", MANIFEST, "--output-path", DEMUX_QZA,
        "--input-format", "PairedEndFastqManifestPhred33V2"])
    if not ok:
        sys.exit("ERROR: import failed.")
    step(Q_DIR / "demux.qzv", "demux summary", [
        "qiime", "demux", "summarize", "--i-data", DEMUX_QZA,
        "--o-visualization", Q_DIR / "demux.qzv"])

    dada2_input = DEMUX_QZA
    if RUN_CUTADAPT:
        if step(TRIMMED_QZA, "cutadapt", [
                "qiime", "cutadapt", "trim-paired",
                "--i-demultiplexed-sequences", DEMUX_QZA,
                "--p-front-f", PRIMER_F, "--p-front-r", PRIMER_R,
                "--p-error-rate", "0.1", "--p-cores", THREADS,
                "--o-trimmed-sequences", TRIMMED_QZA]):
            dada2_input = TRIMMED_QZA
        else:
            print("WARNING: Cutadapt failed — DADA2 will run on untrimmed reads.", file=sys.stderr)
    print(f"   DADA2 input: {dada2_input}")

    ok = step(TABLE_QZA, "DADA2", [
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs", dada2_input,
        "--p-trunc-len-f", TRUNC_LEN_F, "--p-trunc-len-r", TRUNC_LEN_R,
        "--p-trim-left-f", TRIM_LEFT_F, "--p-trim-left-r", TRIM_LEFT_R,
        "--p-max-ee-f", MAX_EE_F, "--p-max-ee-r", MAX_EE_R,
        "--p-n-threads", THREADS,
        "--o-table", TABLE_QZA,
        "--o-representative-sequences", REPSEQS_QZA,
        "--o-denoising-stats", Q_DIR / "denoising-stats.qza"])
    if not ok:
        sys.exit("ERROR: DADA2 failed.")
    step(Q_DIR / "denoising-stats.qzv", "DADA2 stats", [
        "qiime", "metadata", "tabulate",
        "--m-input-file", Q_DIR / "denoising-stats.qza",
        "--o-visualization", Q_DIR / "denoising-stats.qzv"])


def require_asvs():
    for p in (TABLE_QZA, REPSEQS_QZA):
        if not p.exists():
            sys.exit(f"ERROR: not found: {p} — run the 'asv' step first.")


# ── 2. Export ASV sequences ──────────────────────────────────────────────────

def export_sequences():
    if FASTA.exists():
        return
    subprocess.run(["qiime", "tools", "export", "--input-path", str(REPSEQS_QZA),
                    "--output-path", str(FASTA.parent)], check=True)


# ── 3. Remote BLAST, one ASV at a time, most abundant first ──────────────────

def run_blast():
    require_asvs()
    export_sequences()
    seqs = read_fasta(FASTA)
    reads = load_table(set(seqs)).sum(axis=0)
    order = sorted(seqs, key=lambda k: float(reads.get(k, 0)), reverse=True)
    done = read_ids(DONE)
    todo = [k for k in order if k not in done]
    print(f"Total {len(seqs)} | done {len(done)} | remaining {len(todo)}", flush=True)

    cmd = ["blastn", "-query", str(ONE), "-db", BLAST_DB, "-remote",
           "-task", "blastn", "-evalue", EVALUE,
           "-max_target_seqs", str(MAX_TARGET_SEQS), "-max_hsps", str(MAX_HSPS),
           "-outfmt", OUTFMT]

    for n, sid in enumerate(todo, 1):
        ONE.write_text(f">{sid}\n{seqs[sid]}\n")
        print(f"[{n}/{len(todo)}] {sid} reads={int(reads.get(sid, 0))}", flush=True)
        ok = False
        for _ in range(ATTEMPTS):
            try:
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
                if r.returncode == 0 and r.stdout.strip():
                    append_line(RAW, r.stdout)
                    append_line(DONE, sid)
                    print("   hits", len(r.stdout.strip().splitlines()), flush=True)
                    ok = True
                    break
                print("   empty/error:", (r.stderr or r.stdout)[:300], flush=True)
            except subprocess.TimeoutExpired:
                print(f"   timeout ({TIMEOUT_S} s)", flush=True)
            except Exception as e:  # network or process errors: retry, then record
                print("   failed:", e, flush=True)
            time.sleep(PAUSE_RETRY_S)
        if not ok:
            if sid not in read_ids(FAILED):
                append_line(FAILED, sid)
            print("   recorded as failed; continuing (retried on next run)", flush=True)
        time.sleep(PAUSE_BETWEEN_S)

    print("Done:", len(read_ids(DONE)), "| failed:", len(read_ids(FAILED) - read_ids(DONE)))


# ── 4. Best hit per ASV and ASV x sample table ───────────────────────────────

def build_tables():
    if not RAW.exists() or RAW.stat().st_size == 0:
        sys.exit("No BLAST results yet — run the 'blast' step first.")
    require_asvs()

    hits = pd.read_csv(RAW, sep="\t", header=None, names=HIT_COLS)
    hits.to_csv(OUT_DIR / "nt_all_hits.csv", index=False)
    best = (hits.sort_values(["FeatureID", "bitscore"], ascending=[True, False])
                .groupby("FeatureID").first().reset_index())
    best.to_csv(OUT_DIR / "nt_best_hit_per_ASV.csv", index=False)
    print("ASVs with a hit:", len(best))
    print(best["SciName"].value_counts().head(20).to_string())

    from qiime2 import Artifact
    seqs = Artifact.load(str(REPSEQS_QZA)).view(pd.Series)
    seqs.index = seqs.index.astype(str)
    table = load_table(set(seqs.index))

    # Label samples with a metadata column if one is configured
    if "<" not in LABEL_COLUMN:
        meta = pd.read_csv(METADATA, sep="\t", dtype=str)
        if LABEL_COLUMN not in meta.columns:
            sys.exit(f"ERROR: column '{LABEL_COLUMN}' not in metadata.")
        id2lab = dict(zip(meta.iloc[:, 0], meta[LABEL_COLUMN]))
        table.index = [id2lab.get(i, i) for i in table.index]

    samples = list(table.index)
    df = table.T.copy()
    df.index.name = "FeatureID"
    df = df.reset_index()
    df["Sequence"] = df["FeatureID"].map(seqs).astype(str)
    df["Length"] = df["Sequence"].str.len()
    df["TotalReads"] = df[samples].sum(axis=1)
    df["Detected_in"] = df[samples].apply(
        lambda r: ";".join(s for s in samples if r[s] > 0), axis=1)
    df = df.merge(best, on="FeatureID", how="left")
    out = OUT_DIR / "all_ASVs_nt_BLAST_by_sample.csv"
    df.sort_values("TotalReads", ascending=False).to_csv(out, index=False)
    print("Wrote", out)


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage not in ("all", "asv", "blast", "tables"):
        sys.exit("Usage: python 03_its_ncbi_nt_blast.py [all|asv|blast|tables]")
    check_config()
    if stage in ("all", "asv"):
        make_asvs()
    if stage in ("all", "blast"):
        run_blast()
    if stage in ("all", "tables"):
        build_tables()
