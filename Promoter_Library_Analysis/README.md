# Promoter Library Activity Analysis Pipeline

This pipeline performs quantitative analysis of promoter libraries from NGS data.
It integrates paired-end read merging (FLASH2), sequence filtering, barcode matching, read aggregation, and activity score calculation into a single automated workflow.

---

## Directory Structure

```
Promoter_Library_Analysis/
├── promoter_library_analysis.py      # Main analysis script
├── requirements.txt                  # Python dependencies
├── tools/
│   └── FLASH2/                       # FLASH2 (built locally)
├── data/
│   ├── Raw_fastq/                    # Paired-end FASTQ files (.fastq.gz)
│   ├── PLasN5_Bacode_Promoter.csv    # Barcode mapping table (PLasN5)
│   └── PLuxN6_Bacode_Promoter.csv    # Barcode mapping table (PLuxN6)
└── results/                          # Output directory (auto-generated)
```

---

# System Requirements

* **Operating System**: Linux
* **Python**: 3.12.x
* **FLASH2**: compiled from source

---

# Python Dependencies

Only the following Python packages are required:

| Package   | Version |
| --------- | ------- |
| pandas    | 2.2.2   |
| biopython | 1.83    |

Example `requirements.txt`:

```
pandas==2.2.2
biopython==1.83
```

Install with:

```bash
pip install -r requirements.txt
```

---

# Installation

## 1. Clone the Repository

```bash
git clone https://github.com/<your-username>/Promoter_Library_Analysis.git
cd Promoter_Library_Analysis
```

## 2. Create a Virtual Environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

---

# FLASH2 Setup

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

# Input Data

Place the following files in the appropriate directories:

| File                               | Location          |
| ---------------------------------- | ----------------- |
| `*_R1*.fastq.gz`, `*_R2*.fastq.gz` | `data/Raw_fastq/` |
| `PLasN5_Bacode_Promoter.csv`       | `data/`           |
| `PLuxN6_Bacode_Promoter.csv`       | `data/`           |

### FASTQ naming

The plasmid type is inferred from the FASTQ filename.
Include one of the following strings in the sample name:

```
PLasN5
PLuxN6
```

Example:

```
PLasN5_bin1_R1.fastq.gz
PLasN5_bin1_R2.fastq.gz

PLuxN6_bin1_R1.fastq.gz
PLuxN6_bin1_R2.fastq.gz
```

---

# Running the Pipeline

Activate the virtual environment and execute:

```bash
source .venv/bin/activate
python promoter_library_analysis.py
```

---

# Output Files

All results are written to the `results/` directory.

| File / Directory             | Description                             |
| ---------------------------- | --------------------------------------- |
| `merged_output/`             | FLASH2 merged reads                     |
| `filtered_output/`           | Reads passing promoter sequence filters |
| `csv_output/`                | Converted read tables                   |
| `barcode_results/`           | Barcode matching results                |
| `processed_barcode_results/` | Trimmed barcode sequences               |
| `final_results/`             | Final aggregated outputs                |
| `total_counts_*.csv`         | Aggregated read counts per promoter     |
| `activity_summary_*.csv`     | Final activity score summary            |

---

# Pipeline Overview

The analysis consists of the following steps:

1. **Environment check**
2. **FASTQ pair detection**
3. **Read merging with FLASH2**
4. **Promoter sequence filtering**
5. **Barcode matching**
6. **Barcode removal and sequence aggregation**
7. **Activity score calculation**

---

# License

MIT License

---
