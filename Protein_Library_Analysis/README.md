# Protein Library Activity Analysis Pipeline

This pipeline performs quantitative analysis of protein libraries from NGS data.
It integrates paired-end read merging (FLASH2), reference matching, barcode classification, read aggregation, and activity score calculation into a single automated workflow.

---

## Directory Structure

```
Protein_Library_Analysis/
├── protein_library_analysis.py       # Main analysis script
├── requirements.txt                  # Python dependencies
├── tools/
│   └── FLASH2/                       # FLASH2 (built locally)
├── data/
│   ├── Raw_fastq/                    # Paired-end FASTQ files (.fastq.gz)
│   ├── 120poolsLibrary_reference.csv # Reference sequence file
│   ├── input.csv                     # Sample/bin metadata
│   └── FACS_input_data/              # Optional FACS-related data
└── results/                          # Output directory (auto-generated)
```

---

## System Requirements

* **Operating System**: Linux
* **Python**: 3.12.x (tested with 3.12.12)
* **FLASH2**: compiled from source

---

## Python Dependencies

Only the following Python packages are required:

| Package | Version |
| ------- | ------- |
| numpy   | 2.0.2   |
| pandas  | 2.2.2   |

Example `requirements.txt`:

```
numpy==2.0.2
pandas==2.2.2
```

Install with:

```bash
pip install -r requirements.txt
```

---

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/<your-username>/Protein_Library_Analysis.git
cd Protein_Library_Analysis
```

### 2. Create a Virtual Environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

---

## FLASH2 Setup

FLASH2 is required for merging paired-end reads.

```bash
cd tools
git clone https://github.com/dstreett/FLASH2.git
cd FLASH2
make
cd ../..
```

Verify installation:

```bash
./tools/FLASH2/flash2 --version
```

---

## Input Data

Place the following files in the appropriate directories:

| File                               | Location          |
| ---------------------------------- | ----------------- |
| `*_R1*.fastq.gz`, `*_R2*.fastq.gz` | `data/Raw_fastq/` |
| `120poolsLibrary_reference.csv`    | `data/`           |
| `input.csv`                        | `data/`           |

---

## Running the Pipeline

Activate the virtual environment and execute:

```bash
source .venv/bin/activate
python protein_library_analysis.py
```

---

## Output Files

All results are written to the `results/` directory.

| File                               | Description                        |
| ---------------------------------- | ---------------------------------- |
| `merged_output/`                   | FLASH2 merged reads                |
| `*_mapping_results.csv`            | Reference-matched reads            |
| `*_barcode_matched.csv`            | Barcode-matched reads              |
| `*_barcode_unmatched.csv`          | Unmatched barcode reads            |
| `*_barcode_stats.csv`              | Barcode-level summary statistics   |
| `*_detailed_barcode_ref_stats.csv` | Detailed barcode–reference counts  |
| `total_counts_*.csv`               | Aggregated read counts per variant |
| `activity_summary_*.csv`           | Final activity score summary       |

---

## Pipeline Overview

The analysis consists of the following steps:

1. **Environment check**
2. **FASTQ pair detection**
3. **Read merging with FLASH2**
4. **Exact reference sequence matching**
5. **Barcode classification**
6. **Barcode removal and sequence aggregation**
7. **Activity score calculation**

---

## License

MIT License

---
