#!/usr/bin/env python3
"""
Protein Library Analysis Pipeline
==================================
Pipeline for analyzing protein library activity from NGS data.

Usage:
    python protein_library_analysis.py

Directory Structure:
    Protein_Library_Analysis/
    ├── protein_library_analysis.py   # This script
    ├── tools/
    │   └── FLASH2/                   # FLASH2 (Pre-built)
    ├── data/
    │   ├── Raw_fastq/                # R1/R2 FASTQ files (.fastq.gz)
    │   ├── 120poolsLibrary_reference.csv
    │   ├── input.csv
    │   └── FACS_input_data/          # (FACS-related data)
    └── results/                      # Analysis output (Auto-generated)
        ├── merged_output/
        ├── *_mapping_results.csv
        ├── *_barcode_matched.csv
        ├── *_barcode_unmatched.csv
        ├── *_barcode_stats.csv
        ├── *_detailed_barcode_ref_stats.csv
        ├── total_counts_*.csv
        └── activity_summary_*.csv

Preparation:
    1. Build FLASH2 and place it in tools/FLASH2/:
       cd tools && git clone https://github.com/dstreett/FLASH2.git && cd FLASH2 && make
    2. Place Raw FASTQ files in data/Raw_fastq/
    3. Place 120poolsLibrary_reference.csv in data/
    4. Place input.csv in data/
"""

# ============================================
# Library Imports
# ============================================

# Standard Libraries
import os
import subprocess
import glob
import re
import gzip
import math
import traceback
from pathlib import Path
from collections import defaultdict

# Data Analysis
import numpy as np
import pandas as pd

# ============================================
# Path Settings
# ============================================

# Project Root: Directory where this script is located
PROJECT_ROOT = Path(__file__).resolve().parent

# Tools
FLASH2_DIR = PROJECT_ROOT / "tools" / "FLASH2"
FLASH2_BIN = FLASH2_DIR / "flash2"

# Input Data
DATA_DIR = PROJECT_ROOT / "data"
RAW_FASTQ_DIR = DATA_DIR / "Raw_fastq"
REFERENCE_CSV = DATA_DIR / "120poolsLibrary_reference.csv"
INPUT_CSV = DATA_DIR / "input.csv"
FACS_DATA_DIR = DATA_DIR / "FACS_input_data"

# Output
RESULTS_DIR = PROJECT_ROOT / "results"
MERGED_OUTPUT_DIR = RESULTS_DIR / "merged_output"


# ============================================
# Utility Functions
# ============================================

def natural_sort_key(text):
    """
    Converts a string into a key for natural sorting.
    Example: "#10" → [0, 10] (treats numeric parts as integers)
    """
    def atoi(t):
        return int(t) if t.isdigit() else t
    return [atoi(c) for c in re.split(r'(\d+)', str(text))]


def run_command(cmd, description=""):
    """Executes a command and handles errors"""
    try:
        if isinstance(cmd, str):
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        else:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"✗ Error: {description}")
        print(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
        print(f"Error details: {e.stderr}")
        return None


# ============================================
# Step 0: Environment Check
# ============================================

def check_environment():
    """Verify presence of required files and tools"""
    print("=" * 60)
    print("Step 0: Environment Check")
    print("=" * 60)

    errors = []

    # FLASH2
    if not FLASH2_BIN.exists():
        errors.append(
            f"FLASH2 binary not found: {FLASH2_BIN}\n"
            f"  → Please install it using the following commands:\n"
            f"     cd {PROJECT_ROOT / 'tools'} && "
            f"git clone https://github.com/dstreett/FLASH2.git && cd FLASH2 && make"
        )
    else:
        print(f"✓ FLASH2: {FLASH2_BIN}")

    # Raw FASTQ Directory
    if not RAW_FASTQ_DIR.exists():
        errors.append(f"Raw FASTQ directory not found: {RAW_FASTQ_DIR}")
    else:
        fastq_count = len(list(RAW_FASTQ_DIR.glob("*.fastq.gz")))
        print(f"✓ Raw FASTQ directory: {RAW_FASTQ_DIR} ({fastq_count} files)")

    # Reference CSV
    if not REFERENCE_CSV.exists():
        errors.append(f"Reference CSV not found: {REFERENCE_CSV}")
    else:
        print(f"✓ Reference CSV: {REFERENCE_CSV}")

    # Input CSV
    if not INPUT_CSV.exists():
        errors.append(f"input.csv not found: {INPUT_CSV}")
    else:
        print(f"✓ input.csv: {INPUT_CSV}")

    if errors:
        print("\n✗ Found the following errors:")
        for e in errors:
            print(f"  - {e}")
        raise RuntimeError("Environment check failed. Please fix the errors above.")

    # Create output directories
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    MERGED_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"✓ Results output directory: {RESULTS_DIR}")

    print("✓ Environment check complete\n")


# ============================================
# Step 1: Search for FASTQ Pairs
# ============================================

def find_fastq_pairs():
    """Search for R1/R2 pairs in the Raw FASTQ directory"""
    print("=" * 60)
    print("Step 1: Searching for FASTQ files...")
    print("=" * 60)

    r1_files = []
    patterns = [
        str(RAW_FASTQ_DIR / "*_R1_*.fastq.gz"),
        str(RAW_FASTQ_DIR / "*R1*.fastq.gz"),
        str(RAW_FASTQ_DIR / "*_1.fastq.gz"),
    ]
    for pattern in patterns:
        r1_files.extend(glob.glob(pattern))
    r1_files = sorted(list(set(r1_files)))

    pairs = []
    for r1_file in r1_files:
        r2_candidates = [
            r1_file.replace('_R1_', '_R2_'),
            r1_file.replace('R1', 'R2'),
            r1_file.replace('_1.', '_2.'),
        ]
        r2_file = None
        for candidate in r2_candidates:
            if os.path.exists(candidate):
                r2_file = candidate
                break

        if r2_file:
            basename = os.path.basename(r1_file)
            sample_name = re.sub(r'_R1.*|R1.*|_1\..*', '', basename.replace('.fastq.gz', ''))
            pairs.append((sample_name, r1_file, r2_file))
            print(f"✓ Pair confirmed: {sample_name}")

    print(f"Found {len(pairs)} pairs in total.\n")
    return pairs


# ============================================
# Step 2: Merge with FLASH2
# ============================================

def merge_with_flash2(pairs):
    """Merge R1/R2 pairs using FLASH2"""
    if len(pairs) == 0:
        print("Warning: No pairs to merge. Skipping.")
        return []

    print("=" * 60)
    print("Step 2: Merging with FLASH2...")
    print("=" * 60)

    merged_files = []

    for i, (sample_name, r1, r2) in enumerate(pairs, 1):
        print(f"[{i}/{len(pairs)}] Processing: {sample_name}")

        flash_cmd = [
            str(FLASH2_BIN), r1, r2,
            '-o', f"{sample_name}_merged",
            '-d', str(MERGED_OUTPUT_DIR),
            '-M', '330', '-m', '10', '-t', '2', '-z'
        ]

        result = run_command(flash_cmd, f"Merging {sample_name}")

        if result:
            merged_candidates = [
                MERGED_OUTPUT_DIR / f"{sample_name}_merged.extendedFrags.fastq.gz",
                MERGED_OUTPUT_DIR / f"{sample_name}_merged.extendedFrags.fastq",
            ]
            for candidate in merged_candidates:
                if candidate.exists():
                    merged_files.append((sample_name, str(candidate)))
                    print(f"✓ Merge complete: {candidate.name}")
                    break

    print(f"✓ {len(merged_files)} files merged successfully\n")
    return merged_files


# ============================================
# Step 3: Perfect Match Mapping with Reference
# ============================================

def extract_perfect_matches(fastq_file, reference_csv, output_csv):
    """Extract reads that perfectly match reference sequences from FASTQ file"""
    print(f"\nProcessing: {os.path.basename(fastq_file)}")

    # Load Reference CSV
    ref_df = pd.read_csv(reference_csv)
    name_col = next((col for col in ref_df.columns if 'name' in col.lower()), None)
    seq_col = next((col for col in ref_df.columns if 'seq' in col.lower()), None)

    if not name_col or not seq_col:
        print(f"Error: Could not find 'name' or 'seq' columns in CSV")
        print(f"Detected columns: {ref_df.columns.tolist()}")
        return

    # Store reference sequences in a dictionary
    ref_sequences = {}
    for _, row in ref_df.iterrows():
        name = str(row[name_col]).strip()
        seq = str(row[seq_col]).strip().upper()
        ref_sequences[seq] = name

    print(f"✓ Reference sequence count: {len(ref_sequences)}")

    results = []
    total_reads = 0
    matched_reads = 0

    opener = gzip.open if fastq_file.endswith('.gz') else open

    with opener(fastq_file, 'rt') as f:
        while True:
            read_name_line = f.readline()
            if not read_name_line:
                break

            sequence_line = f.readline().strip().upper()
            _ = f.readline()  # +
            _ = f.readline()  # quality

            total_reads += 1
            read_name = read_name_line.strip()[1:].split()[0]

            matched_refs = []
            for ref_seq, ref_name in ref_sequences.items():
                pos = sequence_line.find(ref_seq)
                if pos != -1:
                    match_start = pos + 1  # 1-based
                    match_end = pos + len(ref_seq)
                    matched_refs.append({
                        'ref_name': ref_name,
                        'ref_seq': ref_seq,
                        'start': match_start,
                        'end': match_end,
                    })

            if matched_refs:
                matched_reads += 1
                ref_names = ','.join([m['ref_name'] for m in matched_refs])
                ref_seqs = ','.join([m['ref_seq'] for m in matched_refs])

                first_match = matched_refs[0]
                if len(matched_refs) > 1:
                    position = '; '.join([f"{m['start']}-{m['end']}" for m in matched_refs])
                else:
                    position = f"{first_match['start']}-{first_match['end']}"

                results.append({
                    'ReadName': read_name,
                    'Reference': ref_names,
                    'Position': position,
                    'Sequence': sequence_line,
                    'RefSequence': ref_seqs,
                })

            if total_reads % 10000 == 0:
                print(f"  Processed reads: {total_reads:,} | Matches: {matched_reads:,}", end='\r')

    print(f"\n✓ Processing complete")
    print(f"  Total reads: {total_reads:,}")
    print(f"  Matched reads: {matched_reads:,}")
    if total_reads > 0:
        print(f"  Match rate: {matched_reads / total_reads * 100:.2f}%")

    if results:
        result_df = pd.DataFrame(results)
        result_df.to_csv(output_csv, index=False)
        print(f"✓ Saved results: {output_csv}")
        print(f"  Number of matches saved: {len(result_df):,}")
    else:
        print("Warning: No matching reads found")


def run_mapping(merged_files):
    """Run mapping for all merged files"""
    print("=" * 60)
    print("Step 3: Perfect Match Sequence Extraction")
    print("=" * 60)

    reference_csv = str(REFERENCE_CSV)
    print(f"✓ Reference CSV: {os.path.basename(reference_csv)}")

    # Collect merged FASTQ files
    fastq_files = [f for _, f in merged_files]

    # If merged_files is empty, search in MERGED_OUTPUT_DIR
    if not fastq_files:
        search_patterns = [
            str(MERGED_OUTPUT_DIR / "*_merged.extendedFrags.fastq.gz"),
            str(MERGED_OUTPUT_DIR / "*_merged.extendedFrags.fastq"),
        ]
        for pattern in search_patterns:
            fastq_files.extend(glob.glob(pattern))
        fastq_files = list(set(fastq_files))

    if not fastq_files:
        print("Error: No merged FASTQ files found")
        return

    print(f"✓ Number of FASTQ files found: {len(fastq_files)}")
    for f in fastq_files:
        print(f"  - {f}")

    for fastq_file in fastq_files:
        base_name = os.path.basename(fastq_file)
        base_name = base_name.replace('_merged.extendedFrags.fastq.gz', '')
        base_name = base_name.replace('_merged.extendedFrags.fastq', '')

        output_csv = str(RESULTS_DIR / f"{base_name}_mapping_results.csv")
        extract_perfect_matches(fastq_file, reference_csv, output_csv)

    print("\n✓ Mapping processing complete\n")


# ============================================
# Step 4: Barcode Classification
# ============================================

def run_barcode_classification():
    """Barcode-based read classification"""
    print("=" * 60)
    print("Step 4: Barcode-based Classification Processing")
    print("=" * 60)

    barcode_file = str(INPUT_CSV)
    if not os.path.exists(barcode_file):
        print(f"Error: Barcode file {barcode_file} not found.")
        return

    print(f"Barcode file found: {barcode_file}")

    barcode_df = pd.read_csv(barcode_file)
    print(f"Barcode file columns: {barcode_df.columns.tolist()}")

    # Detect column names
    fw_barcode_candidates = ['Fw_Barcode', 'fw_barcode', 'forward_barcode', 'FwBarcode']
    rv_barcode_candidates = ['Rv_Barcode', 'rv_barcode', 'reverse_barcode', 'RvBarcode']

    fw_col = None
    for candidate in fw_barcode_candidates:
        if candidate in barcode_df.columns:
            fw_col = candidate
            break
    if fw_col is None:
        for col in barcode_df.columns:
            if re.search(r'fw|forward|5|start', col.lower()):
                fw_col = col
                print(f"Estimating forward barcode column: {fw_col}")
                break

    rv_col = None
    for candidate in rv_barcode_candidates:
        if candidate in barcode_df.columns:
            rv_col = candidate
            break
    if rv_col is None:
        for col in barcode_df.columns:
            if re.search(r'rv|reverse|3|end', col.lower()):
                rv_col = col
                print(f"Estimating reverse barcode column: {rv_col}")
                break

    if fw_col is None or rv_col is None:
        print("Error: Barcode columns not found.")
        print(f"Actual columns: {barcode_df.columns.tolist()}")
        return

    print(f"✓ Barcode info loaded: {len(barcode_df)} barcode pairs")
    print(f"  - Forward barcode column: {fw_col}")
    print(f"  - Reverse barcode column: {rv_col}")

    # Load mapping results
    mapping_files = glob.glob(str(RESULTS_DIR / "*_mapping_results.csv"))

    if not mapping_files:
        print(f"Warning: No mapping result files found in {RESULTS_DIR}")
        return

    for mapping_file in mapping_files:
        sample_name = os.path.basename(mapping_file).replace('_mapping_results.csv', '')
        print(f"\nProcessing: {sample_name}")

        try:
            mapping_df = pd.read_csv(mapping_file)
            print(f"Mapping results: {len(mapping_df)} reads")

            barcode_matches = defaultdict(list)
            unmatched = []

            for idx, row in mapping_df.iterrows():
                sequence = row['Sequence']

                if len(sequence) < 16:
                    unmatched.append(row.to_dict())
                    continue

                fw_barcode = sequence[:8]
                rv_barcode = sequence[-8:]

                match_found = False
                for _, bc_row in barcode_df.iterrows():
                    if (fw_barcode == bc_row[fw_col] and
                            rv_barcode == bc_row[rv_col]):
                        barcode_name = bc_row['name']
                        read_data = row.to_dict()
                        read_data['Barcode_Name'] = barcode_name
                        read_data['Fw_Barcode'] = fw_barcode
                        read_data['Rv_Barcode'] = rv_barcode
                        barcode_matches[barcode_name].append(read_data)
                        match_found = True
                        break

                if not match_found:
                    read_data = row.to_dict()
                    read_data['Fw_Barcode'] = fw_barcode
                    read_data['Rv_Barcode'] = rv_barcode
                    unmatched.append(read_data)

            # Output results
            print(f"\nBarcode matching results:")
            sorted_barcodes = sorted(barcode_matches.items(), key=lambda x: x[0])

            all_matches = []
            for barcode_name, reads in sorted_barcodes:
                all_matches.extend(reads)
                print(f"- {barcode_name}: {len(reads)} reads")

            if all_matches:
                matches_df = pd.DataFrame(all_matches)
                matches_output = str(RESULTS_DIR / f"{sample_name}_barcode_matched.csv")
                matches_df.to_csv(matches_output, index=False)
                print(f"✓ Matched reads: {len(all_matches)} ({len(all_matches) / len(mapping_df) * 100:.1f}%)")

            if unmatched:
                unmatched_df = pd.DataFrame(unmatched)
                unmatched_output = str(RESULTS_DIR / f"{sample_name}_barcode_unmatched.csv")
                unmatched_df.to_csv(unmatched_output, index=False)
                print(f"✓ Unmatched reads: {len(unmatched)} ({len(unmatched) / len(mapping_df) * 100:.1f}%)")

            # Barcode statistics
            barcode_stats = []
            for barcode_name, reads in barcode_matches.items():
                if reads:
                    match_df = pd.DataFrame(reads)
                    ref_counts = match_df['Reference'].value_counts().to_dict()
                    top_ref = max(ref_counts.items(), key=lambda x: x[1]) if ref_counts else ('None', 0)
                    barcode_stats.append({
                        'Barcode_Name': barcode_name,
                        'Read_Count': len(reads),
                        'Percent_Total': len(reads) / len(mapping_df) * 100,
                        'Top_Reference': top_ref[0],
                        'Top_Reference_Count': top_ref[1],
                        'Top_Reference_Percent': top_ref[1] / len(reads) * 100 if len(reads) > 0 else 0,
                    })

            if barcode_stats:
                stats_df = pd.DataFrame(barcode_stats)
                stats_df = stats_df.sort_values('Read_Count', ascending=False)
                stats_output = str(RESULTS_DIR / f"{sample_name}_barcode_stats.csv")
                stats_df.to_csv(stats_output, index=False)
                print(f"✓ Barcode stats info: {os.path.basename(stats_output)}")

            # Detailed statistics
            detailed_stats = []
            for barcode_name, reads in barcode_matches.items():
                match_df = pd.DataFrame(reads)
                ref_counts = match_df['Reference'].value_counts().reset_index()
                ref_counts.columns = ['Reference', 'Count']
                for _, ref_row in ref_counts.iterrows():
                    detailed_stats.append({
                        'Barcode_Name': barcode_name,
                        'Reference': ref_row['Reference'],
                        'Read_Count': ref_row['Count'],
                        'Percent_of_Barcode': ref_row['Count'] / len(reads) * 100 if len(reads) > 0 else 0,
                    })

            if detailed_stats:
                detailed_df = pd.DataFrame(detailed_stats)
                detailed_df = detailed_df.sort_values(['Barcode_Name', 'Read_Count'], ascending=[True, False])
                detailed_output = str(RESULTS_DIR / f"{sample_name}_detailed_barcode_ref_stats.csv")
                detailed_df.to_csv(detailed_output, index=False)
                print(f"✓ Detailed barcode-reference stats: {os.path.basename(detailed_output)}")

        except Exception as e:
            print(f"Error: Problem occurred during processing of {mapping_file}.")
            print(f"Details: {str(e)}")
            print(traceback.format_exc())

    print("\n✓ Barcode classification complete\n")


# ============================================
# Step 5: Barcode Removal + Sequence Aggregation + create total_counts
# ============================================

def run_barcode_removal_and_aggregation():
    """Remove barcode sequences, aggregate identical sequences, and create total_counts CSV"""
    print("=" * 60)
    print("Step 5: Barcode Removal and Sequence Aggregation Processing")
    print("=" * 60)

    csv_files = glob.glob(str(RESULTS_DIR / "*_barcode_matched.csv"))

    if not csv_files:
        print(f"Warning: No barcode matched CSV files found in {RESULTS_DIR}")
        return

    print(f"Number of files to process: {len(csv_files)}")

    for csv_file in csv_files:
        print(f"\nProcessing file: {os.path.basename(csv_file)}")

        df = pd.read_csv(csv_file)
        print(f"Number of rows read: {len(df)}")

        df['Original_Sequence'] = df['Sequence']

        # Delete barcode portions
        print("Removing barcode portions from Sequence...")
        modified_sequences = []

        for idx, row in df.iterrows():
            sequence = row['Sequence']
            fw_barcode = row['Fw_Barcode']
            rv_barcode = row['Rv_Barcode']

            if (pd.notna(fw_barcode) and pd.notna(rv_barcode) and
                    isinstance(fw_barcode, str) and isinstance(rv_barcode, str)):
                fw_start = fw_barcode[:8] if len(fw_barcode) >= 8 else fw_barcode
                rv_end = rv_barcode[-8:] if len(rv_barcode) >= 8 else rv_barcode

                if sequence.startswith(fw_start):
                    sequence = sequence[len(fw_start):]
                if sequence.endswith(rv_end):
                    sequence = sequence[:-len(rv_end)]

            modified_sequences.append(sequence)

        df['Sequence'] = modified_sequences

        print(f"Avg sequence length before removal: {df['Original_Sequence'].str.len().mean():.2f}")
        print(f"Avg sequence length after removal: {df['Sequence'].str.len().mean():.2f}")

        df['seq_len'] = df['Sequence'].str.len()

        sequence_counts = df.groupby('Sequence').size().reset_index(name='total')
        barcode_counts = df.groupby(['Sequence', 'Barcode_Name']).size().reset_index(name='count')

        unique_sequences = df.drop_duplicates(subset=['Sequence']).copy()

        barcode_pivot = barcode_counts.pivot(index='Sequence', columns='Barcode_Name', values='count')
        barcode_pivot = barcode_pivot.fillna(0).astype(int)

        merged_df = unique_sequences.merge(sequence_counts, on='Sequence')
        merged_df = merged_df.set_index('Sequence')

        for column in barcode_pivot.columns:
            merged_df[f'Count_{column}'] = barcode_pivot[column]

        merged_df = merged_df.reset_index()

        # Clean columns
        columns_to_drop = ['Original_Sequence', 'Position', 'RefSequence',
                           'Fw_Barcode', 'Rv_Barcode', 'Barcode_Name', 'ReadName']
        existing_cols_to_drop = [col for col in columns_to_drop if col in merged_df.columns]
        df_cleaned = merged_df.drop(columns=existing_cols_to_drop)

        if 'Reference' in df_cleaned.columns:
            df_cleaned['_sort_key'] = df_cleaned['Reference'].apply(natural_sort_key)
            df_cleaned = df_cleaned.sort_values('_sort_key')
            df_cleaned = df_cleaned.drop(columns=['_sort_key'])
            print("✓ Natural sort by Reference complete")

        count_columns = [col for col in df_cleaned.columns if col.startswith('Count_')]
        count_columns_sorted = sorted(count_columns)

        final_columns = ['Reference', 'Sequence', 'seq_len', 'total'] + count_columns_sorted
        final_columns = [col for col in final_columns if col in df_cleaned.columns]
        df_final = df_cleaned[final_columns].reset_index(drop=True)

        df_final.insert(0, 'ReadName', range(1, len(df_final) + 1))

        # Output
        base_name = os.path.basename(csv_file)
        sample_name = base_name.replace('_barcode_matched.csv', '')
        output_name = f'total_counts_{sample_name}.csv'
        output_path = str(RESULTS_DIR / output_name)

        df_final.to_csv(output_path, index=False)

        print(f"\n✓ Processing complete: saved to {output_name}")
        print(f"  Original rows: {len(df):,}")
        print(f"  Aggregated rows: {len(df_final):,}")
        print(f"  Aggregation rate: {(1 - len(df_final) / len(df)) * 100:.2f}%")
        print(f"  Number of Count columns: {len(count_columns_sorted)}")

    # File list
    print("\nGenerated files list:")
    result_files = sorted(glob.glob(str(RESULTS_DIR / "total_counts_*.csv")))
    for result_file in result_files:
        size = os.path.getsize(result_file) / 1024
        print(f"  ✓ {os.path.basename(result_file)} ({size:.1f} KB)")

    print("\n✓ Barcode removal and aggregation complete\n")


# ============================================
# Step 6: Statistics Calculation + create activity_summary
# ============================================

def extract_sample_names_from_input(input_df):
    """Extract sample names from the 'name' column in input.csv"""
    sample_names = set()
    for name in input_df['name']:
        if '_' in name:
            sample_name = name.split('_')[0]
            sample_names.add(sample_name)
    return sorted(list(sample_names))


def calculate_statistics(input_df, total_counts_df, sample_names):
    """Calculate statistics for each read"""
    result_rows = []

    # Calculate total counts for each bin
    bin_total_counts = {}
    for col in total_counts_df.columns:
        if col.startswith('Count_'):
            bin_total_counts[col] = total_counts_df[col].sum()
            print(f"Bin {col} total count: {bin_total_counts[col]}")

    # Identify count columns corresponding to each sample's bins
    sample_bin_columns = {}
    for sample in sample_names:
        sample_bin_columns[sample] = []
        sample_bins = input_df[input_df['name'].str.startswith(f"{sample}_")]['name']
        for bin_name in sample_bins:
            count_col = f"Count_{bin_name}"
            if count_col in total_counts_df.columns:
                sample_bin_columns[sample].append(count_col)
            else:
                print(f"Warning: Count column {count_col} not found in total_counts file")

    for sample, cols in sample_bin_columns.items():
        print(f"Sample {sample} bin count: {len(cols)}")

    # Process each read
    for idx, seq_row in total_counts_df.iterrows():
        res = {
            "ReadName": seq_row["ReadName"],
            "Reference": seq_row["Reference"],
            "Sequence": seq_row["Sequence"],
            "seq_len": seq_row["seq_len"],
            "total": seq_row["total"],
        }

        for sample in sample_names:
            bin_columns = sample_bin_columns[sample]

            sample_total = sum(seq_row[col] for col in bin_columns if col in seq_row)
            res[f"total_{sample}"] = sample_total

            activity_sum = 0.0
            total_weight = 0.0
            weighted_sum_sq = 0.0
            weighted_log_sum = 0.0

            for bin_col in bin_columns:
                if bin_col not in seq_row:
                    continue

                bin_name = bin_col.replace("Count_", "")
                input_row = input_df[input_df['name'] == bin_name]
                if len(input_row) == 0:
                    continue

                gate_min = float(input_row['gate_min'].values[0])
                gate_max = float(input_row['gate_max'].values[0])

                cells_in_gate = 1.0
                if 'cells_in_gate' in input_row.columns:
                    cells_in_gate = float(input_row['cells_in_gate'].values[0])

                rep_value = (gate_min + gate_max) / 2.0
                count = float(seq_row[bin_col])
                total_seq_count = bin_total_counts[bin_col]

                if total_seq_count > 0 and count > 0:
                    weight = (count / total_seq_count) * cells_in_gate
                    activity_sum += rep_value * weight
                    total_weight += weight
                    weighted_sum_sq += weight * (rep_value ** 2)
                    if rep_value > 0:
                        weighted_log_sum += weight * math.log(rep_value)

            if total_weight > 0:
                gauss_mean = activity_sum / total_weight
                gauss_variance = (weighted_sum_sq / total_weight) - (gauss_mean ** 2)
                gauss_variance = max(0.0, gauss_variance)
                geom_mean = math.exp(weighted_log_sum / total_weight) if weighted_log_sum > 0 else 0.0
            else:
                gauss_mean = 0.0
                gauss_variance = 0.0
                geom_mean = 0.0

            res[f"activity_{sample}"] = activity_sum
            res[f"gauss_mean_{sample}"] = gauss_mean
            res[f"gauss_variance_{sample}"] = gauss_variance
            res[f"geom_mean_{sample}"] = geom_mean

        result_rows.append(res)

    return result_rows


def run_activity_summary():
    """Calculate statistics and create activity_summary CSV"""
    print("=" * 60)
    print("Step 6: Statistics Calculation + activity_summary Creation")
    print("=" * 60)

    # Load data
    input_df = pd.read_csv(str(INPUT_CSV))

    total_counts_files = glob.glob(str(RESULTS_DIR / "total_counts_*.csv"))
    if not total_counts_files:
        print("Error: total_counts_*.csv not found")
        return

    total_counts_file = total_counts_files[0]
    total_counts_df = pd.read_csv(total_counts_file)
    total_counts_filename = os.path.basename(total_counts_file)

    print(f"Loaded input.csv with {len(input_df)} rows")
    print(f"Loaded {total_counts_filename} with {len(total_counts_df)} rows")

    # Extract plasmid name
    plasmid_match = re.search(r'total_counts_(.+?)\.csv', total_counts_filename)
    plasmid = plasmid_match.group(1) if plasmid_match else "unknown"
    print(f"Extracted plasmid name: {plasmid}")

    # Extract sample names
    sample_names = extract_sample_names_from_input(input_df)
    print(f"Detected samples from input.csv: {sample_names}")

    # Calculate metrics
    print("\nCalculating statistics...")
    result_rows = calculate_statistics(input_df, total_counts_df, sample_names)

    # Create output headers
    fixed_header = ["ReadName", "Reference", "Sequence", "seq_len", "total"]
    samples_sorted = sorted(sample_names)
    total_headers = [f"total_{s}" for s in samples_sorted]
    activity_headers = [f"activity_{s}" for s in samples_sorted]
    gauss_mean_headers = [f"gauss_mean_{s}" for s in samples_sorted]
    gauss_variance_headers = [f"gauss_variance_{s}" for s in samples_sorted]
    geom_mean_headers = [f"geom_mean_{s}" for s in samples_sorted]

    header = (fixed_header + total_headers + activity_headers +
              gauss_mean_headers + gauss_variance_headers + geom_mean_headers)

    result_df = pd.DataFrame(result_rows)
    valid_headers = [h for h in header if h in result_df.columns]
    result_df = result_df[valid_headers]

    # Remove rows where statistics are 0 or NaN (excluding gauss_variance)
    stat_columns = []
    for sample in sample_names:
        stat_columns.extend([f'gauss_mean_{sample}', f'geom_mean_{sample}'])

    rows_before = len(result_df)

    invalid_rows = result_df[
        (result_df[stat_columns].isna().any(axis=1)) |
        ((result_df[stat_columns] == 0).any(axis=1))
    ].index

    result_df = result_df.drop(invalid_rows)
    rows_after = len(result_df)

    print(f"\nRemoved rows with invalid stats: {rows_before - rows_after} rows deleted ({rows_before} → {rows_after})")

    # Output
    output_filename = f'activity_summary_{plasmid}.csv'
    output_path = str(RESULTS_DIR / output_filename)
    result_df.to_csv(output_path, index=False)

    print(f"\nResults saved to {output_path}")
    print(f"Output file size: {os.path.getsize(output_path) / 1024:.1f}KB")
    print(f"Number of rows: {len(result_df)}")
    print(f"Number of columns: {len(result_df.columns)}")

    # Statistics summary
    for sample in sample_names:
        print(f"\nStats for sample {sample}:")
        for stat_name, col_prefix in [('Gaussian Mean', 'gauss_mean'),
                                       ('Gaussian Variance', 'gauss_variance'),
                                       ('Geometric Mean', 'geom_mean')]:
            col = f'{col_prefix}_{sample}'
            if col in result_df.columns:
                valid_values = result_df[col].dropna()
                if not valid_values.empty:
                    print(f"  {stat_name} range: {valid_values.min():.2f} - {valid_values.max():.2f}")
                    if col_prefix == 'gauss_variance':
                        print(f"  Rows with {stat_name}=0: {(result_df[col] == 0).sum()}")
                    else:
                        print(f"  {stat_name} median: {valid_values.median():.2f}")

    print("\n✓ activity_summary creation complete\n")


# ============================================
# Main Pipeline
# ============================================

def main():
    print("=" * 60)
    print("  Protein Library Analysis Pipeline")
    print("=" * 60)
    print(f"Project Root: {PROJECT_ROOT}\n")

    # Step 0: Environment Check
    check_environment()

    # Step 1: Search for FASTQ Pairs
    pairs = find_fastq_pairs()

    # Step 2: Merge with FLASH2
    merged_files = merge_with_flash2(pairs)

    # Step 3: Perfect Match Mapping with Reference
    run_mapping(merged_files)

    # Step 4: Barcode Classification
    run_barcode_classification()

    # Step 5: Barcode Removal + Sequence Aggregation + create total_counts
    run_barcode_removal_and_aggregation()

    # Step 6: Statistics Calculation + activity_summary Creation
    run_activity_summary()

    # Final result file list
    print("=" * 60)
    print("All processes completed successfully!")
    print("=" * 60)
    print("\nFinal Results Files List:")
    for result_file in sorted(glob.glob(str(RESULTS_DIR / "*.csv"))):
        size = os.path.getsize(result_file) / 1024
        print(f"  ✓ {os.path.basename(result_file)} ({size:.1f} KB)")


if __name__ == "__main__":
    main()