#!/usr/bin/env python3
"""
Promoter Library Analysis Pipeline

A pipeline to analyze promoter library activity from NGS raw data.
It generates activity_summary_*.csv using only reads with 100% barcode matches.

Usage:
    python promoter_library_analysis.py

Directory structure:
    Promoter_Library_Analysis/
    ├── promoter_library_analysis.py
    ├── tools/
    │   └── FLASH2/flash2              # Prebuilt binary
    ├── data/
    │   ├── Raw_fastq/                 # R1/R2 FASTQ files (.fastq.gz / .fastq)
    │   ├── PLasN5_Bacode_Promoter.csv
    │   └── PLuxN6_Bacode_Promoter.csv
    └── results/                       # Automatically generated
        ├── merged_output/
        ├── filtered_output/
        ├── csv_output/
        ├── barcode_results/
        ├── processed_barcode_results/
        └── final_results/

Preparation:
    1. Build FLASH2 and place it at tools/FLASH2/flash2
    2. Place raw FASTQ files (R1/R2 pairs) in data/Raw_fastq/
    3. Place PLasN5_Bacode_Promoter.csv / PLuxN6_Bacode_Promoter.csv in data/
    4. Install dependencies:
         pip install biopython pandas
"""

# ============================================
# Imports
# ============================================
import os
import sys
import subprocess
import glob
import gzip
import csv
import re
import math
import traceback
from pathlib import Path
from collections import defaultdict

import pandas as pd

try:
    from Bio import SeqIO
except ImportError:
    print("Biopython is not installed.")
    print("  pip install biopython")
    sys.exit(1)

# ============================================
# Path settings
# ============================================
BASE_DIR = Path(__file__).resolve().parent
TOOLS_DIR = BASE_DIR / "tools"
FLASH2_PATH = TOOLS_DIR / "FLASH2" / "flash2"
DATA_DIR = BASE_DIR / "data"
RAW_FASTQ_DIR = DATA_DIR / "Raw_fastq"
RESULTS_DIR = BASE_DIR / "results"

MERGED_DIR = RESULTS_DIR / "merged_output"
FILTERED_DIR = RESULTS_DIR / "filtered_output"
CSV_DIR = RESULTS_DIR / "csv_output"
BARCODE_DIR = RESULTS_DIR / "barcode_results"
PROCESSED_DIR = RESULTS_DIR / "processed_barcode_results"
FINAL_DIR = RESULTS_DIR / "final_results"

# Sequence patterns for each plasmid
PLASMID_PATTERNS = {
    "PLasN5": (
        "TTTCTGGAATTCGCGGCCGCTTCTAGAGTTCGAGCCTAGCAAGGGTCCGGGTTCACCGAAACCT"
        "[AGTC]{5}"
        "ATTTGCTAGTTATAAAATTATGAAATTTGCGTAAATTCTTCATACTAGAGGTCGACTGACGACTGGATCCTGTCGGA"
    ),
    "PLuxN6": (
        "TTTCTGGAATTCGCGGCCGCTTCTAGAGTTCGAGCCTAGCAAGGGTCCGGGTTCACACCT"
        "[AGTC]{3}GGATCG[AGTC]{3}"
        "AGGTTTACGCAAGAAAATGGTTTGTTATAGTCGAATAAATACTAGAGGTCGACTGACGACTGGATCCTGTCGGA"
    ),
}


# ============================================
# Utilities
# ============================================
def ensure_dir(path: Path):
    """Create a directory if it does not already exist."""
    path.mkdir(parents=True, exist_ok=True)


def run_command(cmd, description=""):
    """Run a command with error handling."""
    try:
        if isinstance(cmd, str):
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        else:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Error: {description}")
        print(f"    Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
        print(f"    stderr : {e.stderr}")
        return None


def extract_sample_name(filename):
    """
    Extract a standardized sample_name (e.g., CLasR-PLasN5)
    from various filename patterns.

    Supported patterns:
      - CLasR-PLasN5_S381_L001_R1_001.fastq.gz
      - CLasR-PLasN5_merged.extendedFrags.fastq.gz
      - CLasR-PLasN5_S381_L001_filtered.csv
      - CLasR-PLasN5.csv
    """
    basename = os.path.basename(filename)
    # Extract only the {prefix}-{plasmid} part
    m = re.match(r"^([A-Za-z0-9]+-[A-Za-z0-9]+)", basename)
    if m:
        return m.group(1)
    return None


def extract_prefix_and_plasmid(sample_name):
    """
    Split sample_name (e.g., CLasR-PLasN5) into prefix and plasmid.
    """
    if sample_name and "-" in sample_name:
        parts = sample_name.split("-", 1)
        return parts[0], parts[1]
    return None, None


def detect_plasmids():
    """
    Retrieve a list of plasmid names from *_Bacode_Promoter.csv files in data/.
    """
    barcode_files = glob.glob(str(DATA_DIR / "*_Bacode_Promoter.csv"))
    plasmids = []
    for f in barcode_files:
        basename = os.path.basename(f)
        m = re.match(r"^(.+?)_Bacode_Promoter\.csv$", basename)
        if m:
            plasmids.append(m.group(1))
    return sorted(plasmids)


def get_barcode_file(plasmid):
    """Return the path to the barcode file corresponding to a plasmid."""
    return DATA_DIR / f"{plasmid}_Bacode_Promoter.csv"


# ============================================
# Step 1: Merge paired-end reads with FLASH2
# ============================================
def step1_flash2_merge():
    """Detect R1/R2 FASTQ pairs and merge them using FLASH2."""
    print("\n" + "=" * 60)
    print("Step 1: Merge paired-end reads with FLASH2")
    print("=" * 60)

    if not FLASH2_PATH.exists():
        print(f"  ✗ FLASH2 not found: {FLASH2_PATH}")
        print("    Please place the prebuilt binary at tools/FLASH2/flash2.")
        sys.exit(1)

    ensure_dir(MERGED_DIR)

    # Search for R1 files
    all_fastq = glob.glob(str(RAW_FASTQ_DIR / "*.fastq.gz")) + \
                glob.glob(str(RAW_FASTQ_DIR / "*.fastq"))
    r1_files = sorted([f for f in all_fastq if "_R1_" in os.path.basename(f)])

    # Find R1/R2 pairs
    pairs = []
    for r1_file in r1_files:
        r2_file = r1_file.replace("_R1_", "_R2_")
        if os.path.exists(r2_file):
            sample_name = extract_sample_name(r1_file)
            pairs.append((sample_name, r1_file, r2_file))
            print(f"  ✓ Pair detected: {sample_name}")
        else:
            print(f"  ✗ R2 not found: {os.path.basename(r1_file)}")

    print(f"\n  Total pairs found: {len(pairs)}")

    if len(pairs) == 0:
        print("  No input pairs found. Exiting.")
        sys.exit(1)

    # Run FLASH2 merge
    merged_files = []
    for i, (sample_name, r1, r2) in enumerate(pairs, 1):
        print(f"\n  [{i}/{len(pairs)}] Merging: {sample_name}")

        flash_cmd = [
            str(FLASH2_PATH), r1, r2,
            "-o", f"{sample_name}_merged",
            "-d", str(MERGED_DIR),
            "-M", "160", "-m", "10", "-t", "2", "-z"
        ]

        result = run_command(flash_cmd, f"Merging {sample_name}")

        if result:
            # Search for output files
            for candidate in [
                MERGED_DIR / f"{sample_name}_merged.extendedFrags.fastq.gz",
                MERGED_DIR / f"{sample_name}_merged.extendedFrags.fastq",
            ]:
                if candidate.exists():
                    merged_files.append((sample_name, str(candidate)))
                    print(f"  ✓ Merge completed: {candidate.name}")
                    break

    print(f"\n  ✓ Merging completed for {len(merged_files)} files")
    return merged_files


# ============================================
# Step 2: Filter by pattern (with plasmid assignment)
# ============================================
def step2_filter_by_pattern():
    """
    Extract only reads matching the expected sequence pattern
    from merged FASTQ files.

    The plasmid type is inferred from the filename, and only the
    corresponding pattern is used for matching.
    """
    print("\n" + "=" * 60)
    print("Step 2: Sequence pattern matching and filtering")
    print("=" * 60)

    ensure_dir(FILTERED_DIR)

    merged_files = glob.glob(str(MERGED_DIR / "*_merged.extendedFrags.fastq.gz")) + \
                   glob.glob(str(MERGED_DIR / "*_merged.extendedFrags.fastq"))

    print(f"  Input files: {len(merged_files)}\n")

    total_stats = {"total": 0, "matched": 0}
    for pname in PLASMID_PATTERNS:
        total_stats[pname] = 0

    for i, input_file in enumerate(merged_files, 1):
        basename = os.path.basename(input_file)
        # Infer plasmid from sample_name
        sample_name = extract_sample_name(basename)
        _, plasmid = extract_prefix_and_plasmid(sample_name)

        print(f"  [{i}/{len(merged_files)}] {basename}")
        print(f"    Sample name: {sample_name}, inferred plasmid: {plasmid}")

        # Select the relevant pattern
        if plasmid and plasmid in PLASMID_PATTERNS:
            active_patterns = {plasmid: PLASMID_PATTERNS[plasmid]}
        else:
            # If plasmid is unknown, search all patterns
            print(f"    ⚠ Plasmid unknown. Searching all patterns.")
            active_patterns = PLASMID_PATTERNS

        # Output filename: {sample_name}_filtered.fastq.gz
        output_path = FILTERED_DIR / f"{sample_name}_filtered.fastq.gz"

        # Open input file
        if input_file.endswith(".gz"):
            fh_in = gzip.open(input_file, "rt")
        else:
            fh_in = open(input_file, "r")

        fh_out = gzip.open(str(output_path), "wt")

        stats = {"total": 0, "matched": 0}
        for pname in active_patterns:
            stats[pname] = 0

        for record in SeqIO.parse(fh_in, "fastq"):
            stats["total"] += 1
            sequence = str(record.seq)

            matched = False
            for pname, pattern in active_patterns.items():
                if re.search(pattern, sequence):
                    stats[pname] = stats.get(pname, 0) + 1
                    matched = True
                    break

            if matched:
                SeqIO.write(record, fh_out, "fastq")
                stats["matched"] += 1

        fh_in.close()
        fh_out.close()

        # Report
        if stats["total"] > 0:
            pct = stats["matched"] / stats["total"] * 100
        else:
            pct = 0
        print(f"    Total reads: {stats['total']:,}")
        print(f"    Matched    : {stats['matched']:,} ({pct:.2f}%)")
        for pname in active_patterns:
            print(f"      {pname}: {stats.get(pname, 0):,}")

        # Aggregate
        total_stats["total"] += stats["total"]
        total_stats["matched"] += stats["matched"]
        for pname in active_patterns:
            total_stats[pname] = total_stats.get(pname, 0) + stats.get(pname, 0)

    print(f"\n  --- Overall summary ---")
    print(f"  Total reads: {total_stats['total']:,}")
    if total_stats["total"] > 0:
        print(f"  Matched    : {total_stats['matched']:,} "
              f"({total_stats['matched']/total_stats['total']*100:.2f}%)")
    for pname in PLASMID_PATTERNS:
        print(f"    {pname}: {total_stats.get(pname, 0):,}")


# ============================================
# Step 3: Convert FASTQ to CSV
# ============================================
def step3_fastq_to_csv():
    """Convert filtered FASTQ files to CSV format (Sequence_ID, Sequence)."""
    print("\n" + "=" * 60)
    print("Step 3: Convert FASTQ to CSV")
    print("=" * 60)

    ensure_dir(CSV_DIR)

    # Process only _filtered.fastq.gz files generated in Step 2
    filtered_files = sorted(
        glob.glob(str(FILTERED_DIR / "*_filtered.fastq.gz")) +
        glob.glob(str(FILTERED_DIR / "*_filtered.fastq"))
    )

    print(f"  Files to convert: {len(filtered_files)}\n")

    for i, input_file in enumerate(filtered_files, 1):
        basename = os.path.basename(input_file)
        # Extract sample_name (e.g., CLasR-PLasN5_filtered.fastq.gz → CLasR-PLasN5)
        sample_name = basename.split("_filtered")[0]
        output_file = CSV_DIR / f"{sample_name}.csv"

        print(f"  [{i}/{len(filtered_files)}] {basename} → {output_file.name}")

        if input_file.endswith(".gz"):
            opener = gzip.open(input_file, "rt")
        else:
            opener = open(input_file, "r")

        count = 0
        with open(str(output_file), "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Sequence_ID", "Sequence"])

            with opener as f:
                line_count = 0
                seq_id = ""
                sequence = ""
                for line in f:
                    line = line.strip()
                    line_count += 1
                    if line_count == 1:
                        seq_id = line[1:]  # Remove '@'
                    elif line_count == 2:
                        sequence = line
                    elif line_count == 4:
                        writer.writerow([seq_id, sequence])
                        count += 1
                        line_count = 0

        print(f"    ✓ Conversion completed: {count:,} reads")


# ============================================
# Step 4: Barcode matching (100% exact match only)
# ============================================
def step4_barcode_matching():
    """
    Match the first 8 nt and last 8 nt of each read against barcode pairs
    using the plasmid-specific Bacode_Promoter.csv file.
    """
    print("\n" + "=" * 60)
    print("Step 4: Barcode matching (100% exact match)")
    print("=" * 60)

    ensure_dir(BARCODE_DIR)

    # Detect available plasmids
    plasmids = detect_plasmids()
    print(f"  Detected plasmids: {plasmids}")

    # Load barcode information for each plasmid
    barcode_data = {}
    for plasmid in plasmids:
        bc_file = get_barcode_file(plasmid)
        if bc_file.exists():
            df = pd.read_csv(str(bc_file))
            # Normalize column names
            col_map = {}
            for col in df.columns:
                lower = col.lower()
                if re.search(r"fw.*ba[r]?code", lower):
                    col_map[col] = "Fw_Barcode"
                elif re.search(r"rv.*ba[r]?code", lower):
                    col_map[col] = "Rv_Barcode"
            if col_map:
                df = df.rename(columns=col_map)
            barcode_data[plasmid] = df
            print(f"  ✓ {plasmid}: loaded {len(df)} barcode pairs")
            print(f"    Columns: {list(df.columns)}")
        else:
            print(f"  ✗ {bc_file.name} not found")

    # Process CSV files ({sample_name}.csv only, i.e., files containing '-')
    csv_files = sorted(glob.glob(str(CSV_DIR / "*.csv")))
    # Exclude _filtered.csv and _stats.csv; keep only {prefix}-{plasmid}.csv
    csv_files = [f for f in csv_files
                 if not os.path.basename(f).startswith(".")
                 and "_filtered" not in os.path.basename(f)
                 and "_stats" not in os.path.basename(f)
                 and "-" in os.path.basename(f)]

    print(f"\n  CSV files to process: {len(csv_files)}")

    for csv_file in csv_files:
        basename = os.path.basename(csv_file)
        sample_name = basename.replace(".csv", "")
        prefix, plasmid = extract_prefix_and_plasmid(sample_name)

        print(f"\n  Processing: {sample_name} (prefix={prefix}, plasmid={plasmid})")

        if plasmid not in barcode_data:
            print(f"    ✗ No barcode information for plasmid {plasmid}. Skipping.")
            continue

        barcode_df = barcode_data[plasmid]

        # Read input CSV and confirm columns
        df = pd.read_csv(csv_file)
        print(f"    CSV columns: {list(df.columns)}")

        if "Sequence" not in df.columns:
            print(f"    ✗ 'Sequence' column not found. Skipping.")
            print(f"      (Columns: {list(df.columns)})")
            continue

        total_reads = len(df)
        print(f"    Total reads: {total_reads:,}")

        # Build barcode dictionary for faster lookup
        barcode_dict = {}
        for _, bc_row in barcode_df.iterrows():
            key = (str(bc_row["Fw_Barcode"]), str(bc_row["Rv_Barcode"]))
            barcode_dict[key] = bc_row["name"]

        barcode_matches = defaultdict(list)
        unmatched = []

        for _, row in df.iterrows():
            sequence = row["Sequence"]

            if not isinstance(sequence, str) or len(sequence) < 16:
                unmatched.append(row.to_dict())
                continue

            fw_bc = sequence[:8]
            rv_bc = sequence[-8:]

            barcode_name = barcode_dict.get((fw_bc, rv_bc))

            if barcode_name is not None:
                read_data = row.to_dict()
                read_data["Barcode_Name"] = barcode_name
                read_data["Fw_Barcode"] = fw_bc
                read_data["Rv_Barcode"] = rv_bc
                barcode_matches[barcode_name].append(read_data)
            else:
                read_data = row.to_dict()
                read_data["Fw_Barcode"] = fw_bc
                read_data["Rv_Barcode"] = rv_bc
                unmatched.append(read_data)

        # Output results
        all_matches = []
        for barcode_name, reads in sorted(barcode_matches.items()):
            all_matches.extend(reads)
            print(f"      {barcode_name}: {len(reads):,} reads")

        if all_matches:
            matches_df = pd.DataFrame(all_matches)
            out_path = BARCODE_DIR / f"{sample_name}_barcode_matched.csv"
            matches_df.to_csv(str(out_path), index=False)
            pct = len(all_matches) / total_reads * 100 if total_reads > 0 else 0
            print(f"    ✓ Matched  : {len(all_matches):,} ({pct:.1f}%) → {out_path.name}")

        if unmatched:
            unmatched_df = pd.DataFrame(unmatched)
            out_path = BARCODE_DIR / f"{sample_name}_barcode_unmatched.csv"
            unmatched_df.to_csv(str(out_path), index=False)
            pct = len(unmatched) / total_reads * 100 if total_reads > 0 else 0
            print(f"    ✓ Unmatched: {len(unmatched):,} ({pct:.1f}%) → {out_path.name}")


# ============================================
# Step 5: Trim barcodes and merge identical sequences
# ============================================
def step5_trim_and_merge():
    """
    Trim barcode regions (first 8 nt and last 8 nt) and
    merge identical sequences by counting occurrences.
    """
    print("\n" + "=" * 60)
    print("Step 5: Barcode trimming and identical-sequence merging")
    print("=" * 60)

    ensure_dir(PROCESSED_DIR)

    csv_files = sorted(glob.glob(str(BARCODE_DIR / "*_barcode_matched.csv")))

    for csv_file in csv_files:
        basename = os.path.basename(csv_file)
        print(f"\n  Processing: {basename}")

        df = pd.read_csv(csv_file)
        original_len = len(df)

        # Trim barcode regions (remove first 8 nt and last 8 nt)
        trimmed = []
        for _, row in df.iterrows():
            seq = row["Sequence"]
            if isinstance(seq, str) and len(seq) > 16:
                seq = seq[8:-8]
            trimmed.append(seq)

        df["Sequence"] = trimmed

        print(f"    Mean sequence length after trimming: {df['Sequence'].str.len().mean():.1f}")

        # Add sequence length column
        df["seq_len"] = df["Sequence"].str.len()

        # Merge identical sequences
        barcode_counts = df.groupby(["Sequence", "Barcode_Name"]).size().reset_index(name="count")
        barcode_pivot = barcode_counts.pivot(index="Sequence", columns="Barcode_Name", values="count")
        barcode_pivot = barcode_pivot.fillna(0).astype(int)

        # Total counts for each Sequence
        seq_total = df.groupby("Sequence").size().reset_index(name="total")

        # seq_len
        seq_len_df = df.groupby("Sequence")["seq_len"].first().reset_index()

        # Merge
        merged = seq_total.merge(seq_len_df, on="Sequence")
        merged = merged.merge(barcode_pivot, on="Sequence")

        # Reassign ReadName
        merged.insert(0, "ReadName", range(1, len(merged) + 1))

        # Add Count_ prefix
        rename_dict = {}
        for col in barcode_pivot.columns:
            if col in merged.columns:
                rename_dict[col] = f"Count_{col}"
        merged = merged.rename(columns=rename_dict)

        # Reorder columns
        count_cols = sorted([c for c in merged.columns if c.startswith("Count_")])
        final_cols = ["ReadName", "Sequence", "seq_len", "total"] + count_cols
        merged = merged[final_cols]

        # Save
        out_name = basename.replace("_barcode_matched.csv", "_barcode_matched_merged_same_seq.csv")
        out_path = PROCESSED_DIR / out_name
        merged.to_csv(str(out_path), index=False)

        reduction = (1 - len(merged) / original_len) * 100 if original_len > 0 else 0
        print(f"    Original: {original_len:,} rows → Merged: {len(merged):,} rows "
              f"(merge rate: {reduction:.1f}%)")
        print(f"    ✓ Saved: {out_name}")


# ============================================
# Step 6: Combine files by plasmid and generate total_counts
# ============================================
def step6_combine_total_counts():
    """
    Merge processed_barcode_results files by plasmid and
    generate total_counts_{plasmid}.csv.
    """
    print("\n" + "=" * 60)
    print("Step 6: Combine files by plasmid → total_counts")
    print("=" * 60)

    ensure_dir(FINAL_DIR)

    csv_files = sorted(glob.glob(str(PROCESSED_DIR / "*_barcode_matched_merged_same_seq.csv")))

    # Group files by plasmid
    plasmid_files = defaultdict(list)
    for csv_file in csv_files:
        basename = os.path.basename(csv_file)
        sample_name = extract_sample_name(basename)
        if sample_name:
            prefix, plasmid = extract_prefix_and_plasmid(sample_name)
            if prefix and plasmid:
                plasmid_files[plasmid].append((prefix, csv_file))
                print(f"  File detected: {basename} → prefix={prefix}, plasmid={plasmid}")
            else:
                print(f"  ⚠ Failed to parse sample_name: {sample_name}")
        else:
            print(f"  ⚠ Failed to parse filename: {basename}")

    # Merge within each plasmid
    for plasmid, file_list in sorted(plasmid_files.items()):
        print(f"\n  --- Combining {plasmid} ({len(file_list)} files) ---")

        combined_df = None

        for prefix, csv_file in file_list:
            df = pd.read_csv(csv_file)

            # Rename Count_P* columns: Count_P1 → {prefix}-P1
            rename_dict = {}
            for col in df.columns:
                if col.startswith("Count_"):
                    new_name = f"{prefix}-{col[6:]}"  # Remove Count_
                    rename_dict[col] = new_name
            df = df.rename(columns=rename_dict)

            # Drop ReadName (will be reassigned later)
            if "ReadName" in df.columns:
                df = df.drop(columns=["ReadName"])

            if combined_df is None:
                combined_df = df.copy()
            else:
                # Align columns before concatenation
                for col in df.columns:
                    if col not in combined_df.columns:
                        combined_df[col] = 0
                for col in combined_df.columns:
                    if col not in df.columns:
                        df[col] = 0
                combined_df = pd.concat([combined_df, df], ignore_index=True)

        if combined_df is None or len(combined_df) == 0:
            print(f"    ⚠ No data available for {plasmid}")
            continue

        # Merge identical sequences by Sequence
        count_pattern = re.compile(r"^.+-P\d+$")
        count_cols = [c for c in combined_df.columns if count_pattern.match(c)]

        agg_dict = {"seq_len": "first", "total": "sum"}
        for col in count_cols:
            agg_dict[col] = "sum"

        result_df = combined_df.groupby("Sequence", as_index=False).agg(agg_dict)

        # Recalculate total as the sum of count columns
        if count_cols:
            result_df["total"] = result_df[count_cols].sum(axis=1)

        # Reassign ReadName
        result_df.insert(0, "ReadName", range(1, len(result_df) + 1))

        # Reorder columns
        fixed_cols = ["ReadName", "Sequence", "seq_len", "total"]
        other_cols = sorted([c for c in result_df.columns if c not in fixed_cols])
        result_df = result_df[fixed_cols + other_cols]

        # Save
        out_path = FINAL_DIR / f"total_counts_{plasmid}.csv"
        result_df.to_csv(str(out_path), index=False)
        print(f"    ✓ Saved: {out_path.name} ({len(result_df):,} rows, {len(result_df.columns)} columns)")


# ============================================
# Step 7: Calculate activity and generate activity_summary
# ============================================
def step7_calculate_activity():
    """
    Calculate activity statistics for each read using
    total_counts_{plasmid}.csv and {plasmid}_Bacode_Promoter.csv,
    and generate activity_summary_{plasmid}.csv.
    """
    print("\n" + "=" * 60)
    print("Step 7: Calculate activity → activity_summary")
    print("=" * 60)

    plasmids = detect_plasmids()
    print(f"  Target plasmids: {plasmids}")

    for plasmid in plasmids:
        total_counts_file = FINAL_DIR / f"total_counts_{plasmid}.csv"
        barcode_file = get_barcode_file(plasmid)

        if not total_counts_file.exists():
            print(f"\n  ⚠ {total_counts_file.name} not found. Skipping.")
            continue
        if not barcode_file.exists():
            print(f"\n  ⚠ {barcode_file.name} not found. Skipping.")
            continue

        print(f"\n  --- {plasmid} ---")

        total_df = pd.read_csv(str(total_counts_file))
        bacode_df = pd.read_csv(str(barcode_file))

        print(f"    total_counts      : {len(total_df):,} rows")
        print(f"    Bacode_Promoter   : {len(bacode_df)} barcodes")

        # Detect sample names
        sample_names = []
        for col in bacode_df.columns:
            if col.endswith("_cells_in_gate"):
                sample = col.replace("_cells_in_gate", "")
                sample_names.append(sample)
        print(f"    Samples: {sample_names}")

        # Map sample names to total_counts prefixes
        tc_prefixes = set()
        for col in total_df.columns:
            m = re.match(r"^(.+)-P\d+$", col)
            if m:
                tc_prefixes.add(m.group(1))
        print(f"    total_counts prefixes: {sorted(tc_prefixes)}")

        sample_mapping = {}
        for sample in sample_names:
            if sample in tc_prefixes:
                sample_mapping[sample] = sample
            else:
                # Try partial match
                for tc_prefix in tc_prefixes:
                    if tc_prefix in sample or sample in tc_prefix:
                        sample_mapping[sample] = tc_prefix
                        print(f"    Mapping: {sample} → {tc_prefix}")
                        break
                if sample not in sample_mapping:
                    print(f"    ⚠ No matching count columns found for sample {sample}")

        # Precompute total counts for each bin
        bin_total_counts = {}
        for col in total_df.columns:
            if re.match(r".+-P\d+$", col):
                bin_total_counts[col] = total_df[col].sum()

        # Calculate activity
        result_rows = []

        for _, seq_row in total_df.iterrows():
            res = {
                "ReadName": seq_row["ReadName"],
                "Sequence": seq_row["Sequence"],
                "seq_len": seq_row["seq_len"],
                "total": seq_row["total"],
            }

            for sample in sample_names:
                mapped = sample_mapping.get(sample)
                if mapped is None:
                    continue

                bin_columns = [c for c in total_df.columns if c.startswith(f"{mapped}-P")]

                # Total count within the sample
                sample_total = sum(seq_row[c] for c in bin_columns)
                res[f"total_{sample}"] = sample_total

                # Weighted statistics
                activity_sum = 0.0
                total_weight = 0.0
                weighted_sum_sq = 0.0
                weighted_log_sum = 0.0

                for bin_col in bin_columns:
                    bin_name = bin_col.split("-")[1]  # e.g., P1

                    bc_row = bacode_df[bacode_df["name"] == bin_name]
                    if len(bc_row) == 0:
                        continue

                    gate_min = float(bc_row["gate_min"].values[0])
                    gate_max = float(bc_row["gate_max"].values[0])

                    cig_col = f"{sample}_cells_in_gate"
                    if cig_col not in bacode_df.columns:
                        continue
                    cells_in_gate = float(bc_row[cig_col].values[0])

                    rep_value = (gate_min + gate_max) / 2.0
                    count = float(seq_row[bin_col])
                    total_seq_count = bin_total_counts.get(bin_col, 0)

                    if total_seq_count > 0 and count > 0:
                        weight = (count / total_seq_count) * cells_in_gate
                        activity_sum += rep_value * weight
                        total_weight += weight
                        weighted_sum_sq += weight * (rep_value ** 2)
                        if rep_value > 0:
                            weighted_log_sum += weight * math.log(rep_value)

                # Statistics
                if total_weight > 0:
                    gauss_mean = activity_sum / total_weight
                    gauss_variance = (weighted_sum_sq / total_weight) - (gauss_mean ** 2)
                    geom_mean = math.exp(weighted_log_sum / total_weight) if weighted_log_sum != 0 else 0.0
                else:
                    gauss_mean = 0.0
                    gauss_variance = 0.0
                    geom_mean = 0.0

                res[f"activity_{sample}"] = activity_sum
                res[f"gauss_mean_{sample}"] = gauss_mean
                res[f"gauss_variance_{sample}"] = gauss_variance
                res[f"geom_mean_{sample}"] = geom_mean

            result_rows.append(res)

        # Convert to DataFrame and reorder columns
        result_df = pd.DataFrame(result_rows)

        fixed = ["ReadName", "Sequence", "seq_len", "total"]
        samples_sorted = sorted(sample_names)
        total_h = [f"total_{s}" for s in samples_sorted]
        activity_h = [f"activity_{s}" for s in samples_sorted]
        gmean_h = [f"gauss_mean_{s}" for s in samples_sorted]
        gvar_h = [f"gauss_variance_{s}" for s in samples_sorted]
        geom_h = [f"geom_mean_{s}" for s in samples_sorted]

        desired_order = fixed + total_h + activity_h + gmean_h + gvar_h + geom_h
        valid_cols = [c for c in desired_order if c in result_df.columns]
        result_df = result_df[valid_cols]

        # Save
        out_path = FINAL_DIR / f"activity_summary_{plasmid}.csv"
        result_df.to_csv(str(out_path), index=False)
        print(f"    ✓ Saved: {out_path.name} ({len(result_df):,} rows)")


# ============================================
# Main execution
# ============================================
def main():
    print("=" * 60)
    print("  Promoter Library Analysis Pipeline")
    print("=" * 60)
    print(f"  Base directory : {BASE_DIR}")
    print(f"  FLASH2         : {FLASH2_PATH}")
    print(f"  Data directory : {DATA_DIR}")
    print(f"  Output dir     : {RESULTS_DIR}")

    # Check required directories
    if not RAW_FASTQ_DIR.exists():
        print(f"\n  ✗ Raw FASTQ directory not found: {RAW_FASTQ_DIR}")
        sys.exit(1)

    plasmids = detect_plasmids()
    if not plasmids:
        print(f"\n  ✗ Bacode_Promoter.csv files not found.")
        print(f"    Please place files such as PLasN5_Bacode_Promoter.csv in data/.")
        sys.exit(1)
    print(f"  Detected plasmids: {plasmids}")

    # Run pipeline
    step1_flash2_merge()
    step2_filter_by_pattern()
    step3_fastq_to_csv()
    step4_barcode_matching()
    step5_trim_and_merge()
    step6_combine_total_counts()
    step7_calculate_activity()

    print("\n" + "=" * 60)
    print("  All processing completed!")
    print(f"  Results are saved in {FINAL_DIR}")

    # List final output files
    final_files = sorted(glob.glob(str(FINAL_DIR / "*.csv")))
    for f in final_files:
        print(f"    - {os.path.basename(f)}")
    print("=" * 60)


if __name__ == "__main__":
    main()