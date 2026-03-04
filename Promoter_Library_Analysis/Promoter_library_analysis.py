#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Promoter Library Analysis Pipeline
====================================
A pipeline to analyze promoter library activity from NGS Raw data.
Generates activity_summary_*.csv focusing only on reads with a 100% barcode match.

Usage:
    python promoter_library_analysis.py

Directory Structure:
    Promoter_Library_Analysis/
    ├── promoter_library_analysis.py      # This script
    ├── tools/
    │   └── FLASH2/                       # FLASH2 (Pre-built)
    ├── data/
    │   ├── Raw_fastq/                    # R1/R2 FASTQ files (.fastq.gz)
    │   ├── PLasN5_Bacode_Promoter.csv    # Barcode mapping for PLasN5
    │   ├── PLuxN6_Bacode_Promoter.csv    # Barcode mapping for PLuxN6
    │   ├── PLasN5_reference.csv          # (Optional for PLasN5)
    │   └── PLuxN6_reference.csv          # (Optional for PLuxN6)
    └── results/                          # Analysis output (Auto-generated)
        ├── merged_output/                # FLASH2 merge results
        ├── filtered_output/              # Pattern-matched FASTQ
        ├── csv_output/                   # FASTQ to CSV conversion results
        ├── barcode_results/              # Barcode matching results
        │   ├── *_barcode_matched.csv
        │   ├── *_barcode_unmatched.csv
        │   └── *_barcode_stats.csv
        ├── processed_barcode_results/    # Trimming and aggregation results
        │   └── *_barcode_matched_marged_same_seq.csv
        └── final_results/                # Final output
            ├── total_counts_PLasN5.csv
            ├── total_counts_PLuxN6.csv
            ├── activity_summary_PLasN5.csv
            └── activity_summary_PLuxN6.csv

Preparation:
    1. Build FLASH2 and place it in tools/FLASH2/:
       cd tools && git clone https://github.com/dstreett/FLASH2.git && cd FLASH2 && make
    2. Place Raw FASTQ files (R1/R2 pairs) in data/Raw_fastq/
    3. Place PLasN5_Bacode_Promoter.csv and PLuxN6_Bacode_Promoter.csv in data/
    4. Place PLasN5_reference.csv / PLuxN6_reference.csv in data/ if necessary.

Processing Steps:
    Step 1: Paired-end merging via FLASH2
    Step 2: Sequence pattern matching + CSV conversion
    Step 3: Barcode matching (100% match) using plasmid-specific barcode files
    Step 4: Removal of Standard sequences
    Step 5: Barcode trimming + aggregation of identical sequences
    Step 6: Removal of unnecessary columns
    Step 7: Integration of all CSVs + total_count aggregation (by plasmid type)
    Step 8: Statistics calculation + activity_summary.csv creation (by plasmid type)
"""

# ============================================
# Library Imports
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
from Bio import SeqIO


# ============================================
# Set Project Base Directories
# ============================================
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
RAW_FASTQ_DIR = DATA_DIR / "Raw_fastq"
RESULTS_DIR = BASE_DIR / "results"
TOOLS_DIR = BASE_DIR / "tools"

MERGED_DIR = RESULTS_DIR / "merged_output"
FILTERED_DIR = RESULTS_DIR / "filtered_output"
CSV_DIR = RESULTS_DIR / "csv_output"
BARCODE_RESULTS_DIR = RESULTS_DIR / "barcode_results"
PROCESSED_BARCODE_DIR = RESULTS_DIR / "processed_barcode_results"
FINAL_RESULTS_DIR = RESULTS_DIR / "final_results"

FLASH2_PATH = TOOLS_DIR / "FLASH2" / "flash2"

# ============================================
# Pattern Definitions (PLasN5 / PLuxN6)
# ============================================
ALL_SEQUENCE_PATTERNS = {
    'PLasN5': (
        'TTTCTGGAATTCGCGGCCGCTTCTAGAGTTCGAGCCTAGCAAGGGTCCGGGTTCACCGAAACCT'
        '[AGTC]{5}'
        'ATTTGCTAGTTATAAAATTATGAAATTTGCGTAAATTCTTCATACTAGAGGTCGACTGACGACTGGATCCTGTCGGA'
    ),
    'PLuxN6': (
        'TTTCTGGAATTCGCGGCCGCTTCTAGAGTTCGAGCCTAGCAAGGGTCCGGGTTCACACCT'
        '[AGTC]{3}GGATCG[AGTC]{3}'
        'AGGTTTACGCAAGAAAATGGTTTGTTATAGTCGAATAAATACTAGAGGTCGACTGACGACTGGATCCTGTCGGA'
    ),
}


# ============================================
# Utility Functions
# ============================================
def ensure_directories():
    """Create necessary directories"""
    for d in [RESULTS_DIR, MERGED_DIR, FILTERED_DIR, CSV_DIR,
              BARCODE_RESULTS_DIR, PROCESSED_BARCODE_DIR, FINAL_RESULTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)


def run_command(cmd, description=""):
    """Execute command and handle errors"""
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


def detect_plasmid_type(sample_name):
    """Automatically detect plasmid type (PLasN5 or PLuxN6) from sample name"""
    name_upper = sample_name.upper()
    if 'PLASN5' in name_upper:
        return 'PLasN5'
    elif 'PLUXN6' in name_upper:
        return 'PLuxN6'
    else:
        return None


def load_barcode_file(plasmid_type):
    """Load the barcode file corresponding to the plasmid type
    
    Args:
        plasmid_type: 'PLasN5' or 'PLuxN6'
    
    Returns:
        pandas DataFrame or None if file not found
    """
    barcode_filename = f"{plasmid_type}_Bacode_Promoter.csv"
    barcode_path = DATA_DIR / barcode_filename
    
    if not barcode_path.exists():
        print(f"  ✗ Error: Barcode file not found: {barcode_path}")
        return None
    
    df = pd.read_csv(barcode_path)
    print(f"  ✓ Barcode file loaded: {barcode_filename} ({len(df)} rows)")
    return df


def detect_barcode_columns(barcode_df):
    """Detect Fw/Rv barcode column names from the barcode DataFrame
    
    Returns:
        tuple (fw_col, rv_col). (None, None) if not found.
    """
    fw_barcode_candidates = ['Fw_Barcode', 'fw_barcode', 'forward_barcode', 'FwBarcode']
    rv_barcode_candidates = ['Rv_Barcode', 'rv_barcode', 'reverse_barcode', 'RvBarcode']

    fw_col = None
    for candidate in fw_barcode_candidates:
        if candidate in barcode_df.columns:
            fw_col = candidate
            break

    rv_col = None
    for candidate in rv_barcode_candidates:
        if candidate in barcode_df.columns:
            rv_col = candidate
            break

    if fw_col is None:
        for col in barcode_df.columns:
            if re.search(r'fw|forward|5|start', col.lower()):
                fw_col = col
                break

    if rv_col is None:
        for col in barcode_df.columns:
            if re.search(r'rv|reverse|3|end', col.lower()):
                rv_col = col
                break

    return fw_col, rv_col


# ============================================
# Step 1: Paired-end merging via FLASH2
# ============================================
def step1_flash2_merge():
    """Merge paired-end reads using FLASH2"""
    print("\n" + "=" * 60)
    print("Step 1: Paired-end Merging via FLASH2")
    print("=" * 60)

    if not FLASH2_PATH.exists():
        print(f"✗ Error: FLASH2 not found at: {FLASH2_PATH}")
        print("  Please build FLASH2 according to the preparation instructions.")
        sys.exit(1)

    print("\nSearching for FASTQ files...")
    r1_files = []
    search_patterns = [
        str(RAW_FASTQ_DIR / '*_R1_*.fastq.gz'),
        str(RAW_FASTQ_DIR / '*R1*.fastq.gz'),
        str(RAW_FASTQ_DIR / '*_1.fastq.gz'),
    ]
    for pattern in search_patterns:
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
            print(f"  ✓ Pair confirmed: {sample_name}")

    print(f"Total of {len(pairs)} pairs found.")

    if len(pairs) == 0:
        print(f"✗ Error: No FASTQ files found. Please check {RAW_FASTQ_DIR}.")
        sys.exit(1)

    print("\nMerging with FLASH2...")
    merged_files = []

    for i, (sample_name, r1, r2) in enumerate(pairs, 1):
        print(f"[{i}/{len(pairs)}] Processing: {sample_name}")

        flash_cmd = [
            str(FLASH2_PATH), r1, r2,
            '-o', f"{sample_name}_merged",
            '-d', str(MERGED_DIR),
            '-M', '160', '-m', '10', '-t', '2', '-z',
        ]

        result = run_command(flash_cmd, f"Merging {sample_name}")

        if result:
            merged_candidates = [
                str(MERGED_DIR / f"{sample_name}_merged.extendedFrags.fastq.gz"),
                str(MERGED_DIR / f"{sample_name}_merged.extendedFrags.fastq"),
            ]

            for candidate in merged_candidates:
                if os.path.exists(candidate):
                    merged_files.append((sample_name, candidate))
                    print(f"  ✓ Merge complete: {os.path.basename(candidate)}")
                    break

    print(f"\n✓ Merging of {len(merged_files)} files completed")
    return merged_files


# ============================================
# Step 2: Sequence Pattern Matching and Extraction
# ============================================
def step2_filter_and_convert():
    """Filter by sequence pattern matching and convert to CSV"""
    print("\n" + "=" * 60)
    print("Step 2: Sequence Pattern Matching + CSV Conversion")
    print("=" * 60)

    def extract_matching_reads(input_fastq, output_dir, sequence_patterns):
        """Extract reads matching specific patterns from a FASTQ file"""
        os.makedirs(output_dir, exist_ok=True)

        basename = os.path.basename(input_fastq)
        if basename.endswith('.gz'):
            output_filename = basename
        else:
            output_filename = basename + '.gz'

        output_path = os.path.join(output_dir, output_filename)

        if input_fastq.endswith('.gz'):
            file_handle = gzip.open(input_fastq, 'rt')
        else:
            file_handle = open(input_fastq, 'r')

        output_file = gzip.open(output_path, 'wt')

        stats = {'total': 0, 'matched': 0, 'PLasN5': 0, 'PLuxN6': 0}

        for record in SeqIO.parse(file_handle, "fastq"):
            stats['total'] += 1
            sequence = str(record.seq)

            matched = False
            for pattern_name, pattern in sequence_patterns.items():
                if re.search(pattern, sequence):
                    stats[pattern_name] += 1
                    matched = True

            if matched:
                SeqIO.write(record, output_file, "fastq")
                stats['matched'] += 1

        file_handle.close()
        output_file.close()

        print(f"  Total reads: {stats['total']:,}")
        if stats['total'] > 0:
            print(f"  Matched reads: {stats['matched']:,} ({stats['matched']/stats['total']*100:.2f}%)")
        else:
            print(f"  Matched reads: 0")
        for pname in sequence_patterns:
            print(f"    - {pname}: {stats[pname]:,}")

        return output_path, stats

    merged_files = glob.glob(str(MERGED_DIR / '*_merged.extendedFrags.fastq.gz'))
    merged_files.extend(glob.glob(str(MERGED_DIR / '*_merged.extendedFrags.fastq')))

    print(f"\nTarget files: {len(merged_files)}")

    total_stats = {'total': 0, 'matched': 0, 'PLasN5': 0, 'PLuxN6': 0}

    for i, input_file in enumerate(merged_files, 1):
        basename = os.path.basename(input_file)
        sample_name = basename.split('_merged')[0]
        print(f"\n[{i}/{len(merged_files)}] Processing: {basename}")

        plasmid_type = detect_plasmid_type(sample_name)

        if plasmid_type is None:
            print(f"  ⚠ Warning: Could not detect plasmid type from filename: {sample_name}")
            print(f"  → Searching with both patterns (PLasN5, PLuxN6).")
            sequence_patterns = ALL_SEQUENCE_PATTERNS.copy()
        else:
            print(f"  ✓ Auto-detected plasmid type: {plasmid_type}")
            sequence_patterns = {plasmid_type: ALL_SEQUENCE_PATTERNS[plasmid_type]}

        output_path, stats = extract_matching_reads(input_file, str(FILTERED_DIR), sequence_patterns)

        for key in total_stats:
            total_stats[key] += stats[key]

    print("\n" + "-" * 40)
    print("Filtering Completion Summary:")
    print(f"  Total reads: {total_stats['total']:,}")
    if total_stats['total'] > 0:
        print(f"  Matched reads: {total_stats['matched']:,} ({total_stats['matched']/total_stats['total']*100:.2f}%)")
    print(f"    - PLasN5: {total_stats['PLasN5']:,}")
    print(f"    - PLuxN6: {total_stats['PLuxN6']:,}")

    # CSV Conversion
    print("\n" + "-" * 40)
    print("Starting CSV Conversion...")

    filtered_files = glob.glob(str(FILTERED_DIR / '*_merged.extendedFrags.fastq.gz'))
    filtered_files.extend(glob.glob(str(FILTERED_DIR / '*_merged.extendedFrags.fastq')))

    for i, input_file in enumerate(filtered_files, 1):
        basename = os.path.basename(input_file)
        sample_name = basename.split('_merged')[0]
        output_file = str(CSV_DIR / f"{sample_name}.csv")

        print(f"[{i}/{len(filtered_files)}] Converting: {basename} → {os.path.basename(output_file)}")

        with open(output_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(['Sequence_ID', 'Sequence'])

            if input_file.endswith('.gz'):
                opener = gzip.open(input_file, 'rt')
            else:
                opener = open(input_file, 'r')

            with opener as f:
                line_count = 0
                seq_id = ''
                sequence = ''

                for line in f:
                    line = line.strip()
                    line_count += 1

                    if line_count == 1:
                        seq_id = line[1:]
                    elif line_count == 2:
                        sequence = line
                    elif line_count == 4:
                        csvwriter.writerow([seq_id, sequence])
                        line_count = 0

        print(f"  ✓ Conversion complete: {output_file}")

    print(f"\n✓ CSV conversion of all {len(filtered_files)} files completed.")


# ============================================
# Step 3: Barcode Matching (100% match only)
#         ★ Uses plasmid-specific barcode files
# ============================================
def step3_barcode_matching():
    """Classification processing based on barcodes (uses specific barcode files per plasmid type)"""
    print("\n" + "=" * 60)
    print("Step 3: Barcode Matching (100% Match)")
    print("=" * 60)

    # Pre-load barcode files by plasmid type
    barcode_data = {}  # {plasmid_type: (barcode_df, fw_col, rv_col, has_reference)}

    for plasmid_type in ['PLasN5', 'PLuxN6']:
        barcode_df = load_barcode_file(plasmid_type)
        if barcode_df is None:
            print(f"  ⚠ Warning: Barcode file for {plasmid_type} not found. "
                  f"{plasmid_type} samples will be skipped.")
            continue

        print(f"  Columns in {plasmid_type}_Bacode_Promoter.csv: {barcode_df.columns.tolist()}")

        fw_col, rv_col = detect_barcode_columns(barcode_df)

        if fw_col is None or rv_col is None:
            print(f"  ✗ Error: Barcode columns for {plasmid_type} not found.")
            print(f"    Actual columns: {barcode_df.columns.tolist()}")
            continue

        has_reference = 'Reference' in barcode_df.columns
        if has_reference:
            print(f"  ✓ {plasmid_type}: 'Reference' column detected")
        else:
            print(f"  ⚠ {plasmid_type}: No 'Reference' column (will remain blank)")

        print(f"  ✓ {plasmid_type}: {len(barcode_df)} barcode pairs")
        print(f"    Forward barcode column: {fw_col}")
        print(f"    Reverse barcode column: {rv_col}")

        barcode_data[plasmid_type] = (barcode_df, fw_col, rv_col, has_reference)

    if not barcode_data:
        print("✗ Error: No valid barcode files found.")
        sys.exit(1)

    # Process CSV files
    csv_files = glob.glob(str(CSV_DIR / '*.csv'))
    csv_files = [f for f in csv_files if not f.endswith('_stats.csv')]

    if not csv_files:
        print("✗ Error: No target CSV files found.")
        return

    print(f"\nNumber of files to process: {len(csv_files)}")

    for csv_file in csv_files:
        sample_name = os.path.basename(csv_file).replace('.csv', '')
        print(f"\nProcessing: {sample_name}")

        # ★ Determine plasmid type from sample name and select corresponding barcodes
        plasmid_type = detect_plasmid_type(sample_name)

        if plasmid_type is None:
            print(f"  ⚠ Warning: Could not detect plasmid type from sample name: {sample_name}")
            print(f"  → Skipping.")
            continue

        if plasmid_type not in barcode_data:
            print(f"  ⚠ Warning: Barcode file for {plasmid_type} not loaded. Skipping.")
            continue

        barcode_df, fw_col, rv_col, has_reference = barcode_data[plasmid_type]
        print(f"  ✓ Using barcode file: {plasmid_type}_Bacode_Promoter.csv ({len(barcode_df)} entries)")

        try:
            df = pd.read_csv(csv_file)
            total_reads = len(df)
            print(f"  Total reads: {total_reads}")

            barcode_matches = defaultdict(list)
            unmatched = []

            for idx, row in df.iterrows():
                sequence = row['Sequence']

                if len(sequence) < 16:
                    unmatched.append(row.to_dict())
                    continue

                fw_barcode = sequence[:8]
                rv_barcode = sequence[-8:]

                match_found = False
                for bc_idx, bc_row in barcode_df.iterrows():
                    if (fw_barcode == bc_row[fw_col] and
                            rv_barcode == bc_row[rv_col]):
                        barcode_name = bc_row['name'] if 'name' in bc_row else f"Barcode_{bc_idx}"
                        reference_name = bc_row['Reference'] if has_reference else ''
                        read_data = row.to_dict()
                        read_data['Barcode_Name'] = barcode_name
                        read_data['Reference'] = reference_name
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
            print(f"\n  Barcode matching results:")
            all_matches = []
            for barcode_name, reads in barcode_matches.items():
                all_matches.extend(reads)
                print(f"    - {barcode_name}: {len(reads)} reads")

            if all_matches:
                matches_df = pd.DataFrame(all_matches)
                matches_output = str(BARCODE_RESULTS_DIR / f"{sample_name}_barcode_matched.csv")
                matches_df.to_csv(matches_output, index=False)
                print(f"  ✓ Matched reads: {len(all_matches)} ({len(all_matches)/total_reads*100:.1f}%)")

            if unmatched:
                unmatched_df = pd.DataFrame(unmatched)
                unmatched_output = str(BARCODE_RESULTS_DIR / f"{sample_name}_barcode_unmatched.csv")
                unmatched_df.to_csv(unmatched_output, index=False)
                print(f"  ✓ Unmatched reads: {len(unmatched)} ({len(unmatched)/total_reads*100:.1f}%)")

            # Statistics per barcode
            barcode_stats = []
            for barcode_name, reads in barcode_matches.items():
                if reads:
                    barcode_stats.append({
                        'Barcode_Name': barcode_name,
                        'Read_Count': len(reads),
                        'Percent_Total': len(reads) / total_reads * 100,
                    })

            if barcode_stats:
                stats_df = pd.DataFrame(barcode_stats)
                stats_df = stats_df.sort_values('Read_Count', ascending=False)
                stats_output = str(BARCODE_RESULTS_DIR / f"{sample_name}_barcode_stats.csv")
                stats_df.to_csv(stats_output, index=False)
                print(f"  ✓ Barcode statistics: {os.path.basename(stats_output)}")

        except Exception as e:
            print(f"  ✗ Error: Problem occurred during processing of {csv_file}.")
            print(f"  Details: {str(e)}")
            print(traceback.format_exc())

    print(f"\n✓ Barcode classification processing completed.")


# ============================================
# Step 4: Removal of Standard Sequences
# ============================================
def step4_remove_standard():
    """Remove Standard sequences"""
    print("\n" + "=" * 60)
    print("Step 4: Removal of Standard Sequences")
    print("=" * 60)

    files = glob.glob(str(BARCODE_RESULTS_DIR / "*_barcode_matched.csv"))

    if not files:
        print("✗ No target files found.")
        return

    print(f"Targeting: {len(files)} files")

    total_removed = 0
    total_remaining = 0

    for file_path in files:
        try:
            df = pd.read_csv(file_path)
            original_rows = len(df)

            if 'Reference' in df.columns:
                df = df[df['Reference'] != 'Standard']
                removed_rows = original_rows - len(df)
            else:
                removed_rows = 0
                print(f"  ⚠ {os.path.basename(file_path)}: Skipping Standard removal because 'Reference' column is missing")

            df.to_csv(file_path, index=False)

            total_removed += removed_rows
            total_remaining += len(df)

            print(f"  {os.path.basename(file_path)}: {removed_rows} rows removed, {len(df)} rows remaining")
        except Exception as e:
            print(f"  ✗ Error: {e}")

    print(f"\n✓ Total of {total_removed} Standard sequences removed. Remaining: {total_remaining} rows")


# ============================================
# Step 5: Barcode Trimming + Identical Sequence Aggregation
# ============================================
def step5_trim_and_merge():
    """Barcode sequence trimming + aggregation of identical base sequences + per-bin counts + length calculation"""
    print("\n" + "=" * 60)
    print("Step 5: Barcode Trimming + Identical Sequence Aggregation")
    print("=" * 60)

    csv_files = glob.glob(str(BARCODE_RESULTS_DIR / '*_barcode_matched.csv'))

    for csv_file in csv_files:
        print(f"\nProcessing: {os.path.basename(csv_file)}")

        df = pd.read_csv(csv_file)
        df['Original_Sequence'] = df['Sequence']

        modified_sequences = []
        for idx, row in df.iterrows():
            sequence = row['Sequence']
            fw_barcode = row['Fw_Barcode']
            rv_barcode = row['Rv_Barcode']

            if (pd.notna(fw_barcode) and pd.notna(rv_barcode)
                    and isinstance(fw_barcode, str) and isinstance(rv_barcode, str)):
                fw_start = fw_barcode[:8] if len(fw_barcode) >= 8 else fw_barcode
                rv_end = rv_barcode[-8:] if len(rv_barcode) >= 8 else rv_barcode

                if sequence.startswith(fw_start):
                    sequence = sequence[len(fw_start):]
                if sequence.endswith(rv_end):
                    sequence = sequence[:-len(rv_end)]

            modified_sequences.append(sequence)

        df['Sequence'] = modified_sequences

        print(f"  Avg Sequence length before barcode removal: {df['Original_Sequence'].str.len().mean():.2f}")
        print(f"  Avg Sequence length after barcode removal: {df['Sequence'].str.len().mean():.2f}")

        df['seq_len'] = df['Sequence'].str.len()

        sequence_counts = df.groupby('Sequence').size().reset_index(name='total')
        barcode_counts = df.groupby(['Sequence', 'Barcode_Name']).size().reset_index(name='count')

        if 'Reference' in df.columns:
            ref_info = df.drop_duplicates(subset=['Sequence'])[['Sequence', 'Reference']].copy()
        else:
            ref_info = df.drop_duplicates(subset=['Sequence'])[['Sequence']].copy()
            ref_info['Reference'] = ''

        barcode_pivot = barcode_counts.pivot(index='Sequence', columns='Barcode_Name', values='count')
        barcode_pivot = barcode_pivot.fillna(0).astype(int)

        seq_len_info = df.drop_duplicates(subset=['Sequence'])[['Sequence', 'seq_len']].copy()

        merged_df = ref_info.merge(sequence_counts, on='Sequence')
        merged_df = merged_df.merge(seq_len_info, on='Sequence')
        merged_df = merged_df.set_index('Sequence')

        for column in barcode_pivot.columns:
            merged_df[f'Count_{column}'] = barcode_pivot[column]

        merged_df = merged_df.reset_index()

        base_cols = ['Sequence', 'Reference', 'seq_len', 'total']
        count_cols = [col for col in merged_df.columns if col.startswith('Count_')]
        other_cols = [col for col in merged_df.columns
                      if col not in base_cols and col not in count_cols]

        new_cols = base_cols + count_cols + other_cols
        new_cols = [c for c in new_cols if c in merged_df.columns]
        final_df = merged_df[new_cols]

        base_name = os.path.basename(csv_file)
        output_name = base_name.replace('_barcode_matched.csv', '_barcode_matched_marged_same_seq.csv')
        output_path = str(PROCESSED_BARCODE_DIR / output_name)

        final_df.to_csv(output_path, index=False)
        print(f"  ✓ Saved: {output_name}")
        print(f"  Original rows: {len(df)}, Aggregated rows: {len(final_df)}")
        if len(df) > 0:
            print(f"  Aggregation rate: {(1 - len(final_df) / len(df)) * 100:.2f}%")


# ============================================
# Step 6: Removal of Unnecessary Columns
# ============================================
def step6_clean_columns():
    """Remove unnecessary columns"""
    print("\n" + "=" * 60)
    print("Step 6: Cleaning Columns")
    print("=" * 60)

    processed_files = glob.glob(str(PROCESSED_BARCODE_DIR / '*_barcode_matched_marged_same_seq.csv'))

    for file_path in processed_files:
        print(f"\nProcessing: {os.path.basename(file_path)}")

        df = pd.read_csv(file_path)

        columns_to_drop = ['Original_Sequence', 'Fw_Barcode', 'Rv_Barcode', 'Barcode_Name', 'Sequence_ID']
        columns_to_drop = [col for col in columns_to_drop if col in df.columns]

        if columns_to_drop:
            df = df.drop(columns=columns_to_drop)
            print(f"  Removed columns: {columns_to_drop}")

        df.insert(0, 'ReadName', range(1, len(df) + 1))
        df.to_csv(file_path, index=False)
        print(f"  ✓ Save complete (Added 'ReadName' column)")

    print(f"\n✓ Column organization for all files completed.")


# ============================================
# Step 7: Integrate all CSVs and Aggregate total_count (by Plasmid Type)
# ============================================
def step7_integrate_counts():
    """Integrate all CSVs by plasmid type and aggregate total_count"""
    print("\n" + "=" * 60)
    print("Step 7: Global CSV Integration + total_count Aggregation (by Plasmid Type)")
    print("=" * 60)

    csv_files = glob.glob(str(PROCESSED_BARCODE_DIR / '*_barcode_matched_marged_same_seq.csv'))

    if not csv_files:
        print("✗ No target files found.")
        return

    pattern = r'(\w+)-([A-Za-z0-9]+)_'

    plasmid_groups = defaultdict(list)
    for csv_file in csv_files:
        basename = os.path.basename(csv_file)
        match = re.search(pattern, basename)

        if match:
            prefix = match.group(1)
            plasmid = match.group(2)
            plasmid_groups[plasmid].append((prefix, csv_file))
            print(f"  {basename} → Plasmid: {plasmid}, Prefix: {prefix}")
        else:
            print(f"  ✗ Could not extract prefix/plasmid from filename: {basename}")

    if not plasmid_groups:
        print("✗ No classifiable files found.")
        return

    print(f"\nDetected plasmid types: {list(plasmid_groups.keys())}")

    for plasmid, file_list in plasmid_groups.items():
        print(f"\n--- Plasmid: {plasmid} ({len(file_list)} files) ---")

        combined_df = None

        for prefix, csv_file in file_list:
            print(f"  Reading: {os.path.basename(csv_file)}")

            df = pd.read_csv(csv_file)

            count_cols = [col for col in df.columns if col.startswith('Count_P')]
            rename_dict = {}
            for col in count_cols:
                new_col_name = f"{prefix}-{col[6:]}"
                rename_dict[col] = new_col_name

            df = df.rename(columns=rename_dict)

            if combined_df is None:
                combined_df = df.copy()
            else:
                for col in rename_dict.values():
                    if col not in combined_df.columns:
                        combined_df[col] = 0

                for col in combined_df.columns:
                    if col not in df.columns and re.match(r'\w+-P\d+', col):
                        df[col] = 0

                combined_df = pd.concat([combined_df, df], ignore_index=True)

        if combined_df is not None and len(combined_df) > 0:
            agg_dict = {}
            if 'Reference' in combined_df.columns:
                agg_dict['Reference'] = 'first'
            agg_dict['seq_len'] = 'first'
            agg_dict['total'] = 'sum'

            for col in combined_df.columns:
                if re.match(r'\w+-P\d+', col):
                    agg_dict[col] = 'sum'

            result_df = combined_df.groupby('Sequence', as_index=False).agg(agg_dict)
            result_df['ReadName'] = range(1, len(result_df) + 1)

            cols = result_df.columns.tolist()
            cols.remove('ReadName')
            result_df = result_df[['ReadName'] + cols]

            output_name = f"total_counts_{plasmid}.csv"
            output_path = str(FINAL_RESULTS_DIR / output_name)
            result_df.to_csv(output_path, index=False)

            print(f"\n  ✓ Integration complete: {output_name}")
            print(f"    Rows: {len(result_df)}, Columns: {len(result_df.columns)}")
        else:
            print(f"\n  ✗ {plasmid}: No data found for integration.")

    print(f"\n✓ Integration for all plasmid types completed.")


# ============================================
# Step 8: Statistics Calculation + activity_summary.csv (by Plasmid Type)
#         ★ Uses plasmid-specific barcode files
# ============================================
def step8_activity_summary():
    """Calculate statistics and create activity_summary.csv per plasmid type"""
    print("\n" + "=" * 60)
    print("Step 8: Statistics Calculation + activity_summary.csv Creation (by Plasmid Type)")
    print("=" * 60)

    # Get all total_counts files
    total_counts_files = glob.glob(str(FINAL_RESULTS_DIR / 'total_counts_*.csv'))
    if not total_counts_files:
        print("✗ Error: total_counts_*.csv not found.")
        sys.exit(1)

    print(f"Detected total_counts files: {len(total_counts_files)}")
    for f in total_counts_files:
        print(f"  - {os.path.basename(f)}")

    # Create activity_summary for each total_counts file
    for total_counts_file in total_counts_files:
        total_counts_filename = os.path.basename(total_counts_file)
        total_counts_df = pd.read_csv(total_counts_file)

        # Extract plasmid name
        plasmid_match = re.search(r'total_counts_(.+?)\.csv', total_counts_filename)
        plasmid = plasmid_match.group(1) if plasmid_match else "unknown"

        print(f"\n{'=' * 40}")
        print(f"Processing: {total_counts_filename} (Plasmid: {plasmid})")
        print(f"  Rows: {len(total_counts_df)}")

        # ★ Load barcode file corresponding to plasmid type
        plasmid_type = detect_plasmid_type(plasmid)
        if plasmid_type is None:
            print(f"  ⚠ Warning: Could not detect plasmid type for: {plasmid}. Skipping.")
            continue

        bacode_df = load_barcode_file(plasmid_type)
        if bacode_df is None:
            print(f"  ✗ Error: {plasmid_type}_Bacode_Promoter.csv not found. Skipping.")
            continue

        print(f"  ✓ Using barcode file: {plasmid_type}_Bacode_Promoter.csv ({len(bacode_df)} rows)")

        # Extract sample names (from _cells_in_gate columns)
        sample_names = []
        for col in bacode_df.columns:
            if '_cells_in_gate' in col:
                sample = col.replace('_cells_in_gate', '')
                sample_names.append(sample)

        print(f"  Detected samples: {sample_names}")

        # Sample mapping (adjust as needed)
        sample_mapping = {'pET16b': 'pET16'}
        print(f"  Sample mapping: {sample_mapping}")

        # Calculate total counts for each bin
        bin_total_counts = {}
        for sample in sample_names:
            mapped_sample = sample_mapping.get(sample, sample)
            bin_columns = [col for col in total_counts_df.columns if col.startswith(f"{mapped_sample}-P")]
            for col in bin_columns:
                bin_total_counts[col] = total_counts_df[col].sum()

        if not bin_total_counts:
            print(f"  ⚠ Warning: No bin columns found corresponding to {plasmid}. Skipping.")
            print(f"    Columns in total_counts: {total_counts_df.columns.tolist()}")
            print(f"    Search patterns: {[f'{sample_mapping.get(s, s)}-P' for s in sample_names]}")
            continue

        # Calculate metrics
        print(f"  Calculating statistics...")
        result_rows = []

        for _, seq_row in total_counts_df.iterrows():
            res = {
                "ReadName": seq_row["ReadName"],
                "Sequence": seq_row["Sequence"],
                "seq_len": seq_row["seq_len"],
                "total": seq_row["total"],
            }

            if 'Reference' in seq_row.index:
                res["Reference"] = seq_row["Reference"]

            for sample in sample_names:
                mapped_sample = sample_mapping.get(sample, sample)
                bin_columns = [col for col in total_counts_df.columns if col.startswith(f"{mapped_sample}-P")]

                sample_total = sum(seq_row[col] for col in bin_columns)
                res[f"total_{sample}"] = sample_total

                activity_sum = 0.0
                total_weight = 0.0
                weighted_sum_sq = 0.0
                weighted_log_sum = 0.0

                for bin_col in bin_columns:
                    bin_name = bin_col.split('-')[1]
                    bacode_row = bacode_df[bacode_df['name'] == bin_name]
                    if len(bacode_row) == 0:
                        continue

                    gate_min = float(bacode_row['gate_min'].values[0])
                    gate_max = float(bacode_row['gate_max'].values[0])

                    cells_in_gate_col = f"{sample}_cells_in_gate"
                    if cells_in_gate_col not in bacode_row.columns:
                        continue
                    cells_in_gate = float(bacode_row[cells_in_gate_col].values[0])

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

        # Prepare output header
        fixed_header = ["ReadName"]
        if 'Reference' in total_counts_df.columns:
            fixed_header.append("Reference")
        fixed_header.extend(["Sequence", "seq_len", "total"])

        samples_sorted = sorted(sample_names)
        total_headers = [f"total_{s}" for s in samples_sorted]
        activity_headers = [f"activity_{s}" for s in samples_sorted]
        gauss_mean_headers = [f"gauss_mean_{s}" for s in samples_sorted]
        gauss_variance_headers = [f"gauss_variance_{s}" for s in samples_sorted]
        geom_mean_headers = [f"geom_mean_{s}" for s in samples_sorted]

        header = (fixed_header + total_headers + activity_headers
                  + gauss_mean_headers + gauss_variance_headers + geom_mean_headers)

        result_df = pd.DataFrame(result_rows)
        valid_headers = [h for h in header if h in result_df.columns]
        result_df = result_df[valid_headers]

        output_filename = f'activity_summary_{plasmid}.csv'
        output_path = str(FINAL_RESULTS_DIR / output_filename)
        result_df.to_csv(output_path, index=False)
        print(f"\n  ✓ Result saved: {output_filename}")
        print(f"    Rows: {len(result_df)}, Columns: {len(result_df.columns)}")

    print(f"\n✓ activity_summary for all plasmid types completed.")


# ============================================
# Main Pipeline
# ============================================
def main():
    """Entry point for executing the entire pipeline"""
    print("=" * 60)
    print(" Promoter Library Analysis Pipeline")
    print("=" * 60)
    print(f"Project Directory: {BASE_DIR}")
    print(f"Data Directory:    {DATA_DIR}")
    print(f"Results Directory: {RESULTS_DIR}")
    print("=" * 60)

    ensure_directories()

    # ★ Check for the presence of plasmid-specific barcode files
    barcode_files_found = []
    for plasmid_type in ['PLasN5', 'PLuxN6']:
        bf = DATA_DIR / f'{plasmid_type}_Bacode_Promoter.csv'
        if bf.exists():
            barcode_files_found.append(bf)
            print(f"  ✓ {bf.name}")
        else:
            print(f"  ⚠ {bf.name} not found")

    if not barcode_files_found:
        print(f"\n✗ Error: Neither PLasN5_Bacode_Promoter.csv nor PLuxN6_Bacode_Promoter.csv "
              f"was found in data/.")
        sys.exit(1)

    raw_fastq_files = list(RAW_FASTQ_DIR.glob('*.fastq.gz'))
    if not raw_fastq_files:
        print(f"\n✗ Error: No FASTQ files found in {RAW_FASTQ_DIR}.")
        sys.exit(1)

    print(f"\nInput File Verification:")
    print(f"  ✓ Barcode files: {len(barcode_files_found)}")
    print(f"  ✓ FASTQ files: {len(raw_fastq_files)}")

    print(f"\nPlasmid Type Pre-verification:")
    for fq in sorted(raw_fastq_files):
        detected = detect_plasmid_type(fq.name)
        status = detected if detected else "Undetectable (will search with both patterns)"
        print(f"  {fq.name} → {status}")

    # Execute Pipeline
    step1_flash2_merge()
    step2_filter_and_convert()
    step3_barcode_matching()
    step4_remove_standard()
    step5_trim_and_merge()
    step6_clean_columns()
    step7_integrate_counts()
    step8_activity_summary()

    print("\n" + "=" * 60)
    print(" Global Processing Completed!")
    print("=" * 60)
    print(f"\nFinal results are saved in the following directory:")
    print(f"  {FINAL_RESULTS_DIR}")

    final_files = list(FINAL_RESULTS_DIR.glob('*.csv'))
    if final_files:
        for f in sorted(final_files):
            print(f"  - {f.name}")
    else:
        print("  (No files found)")


if __name__ == "__main__":
    main()