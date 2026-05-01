# Stringology-Based Motif Discovery for EEG Analysis

This repository contains the implementation of the OPM (Order Preserving Matching) and CTM (Cartesian Tree Matching) based motif discovery framework described in:

> Dahan, A., & Ghazawi, S. (2026). *Stringology-Based Motif Discovery for Electrophysiological Time Series: A Framework for Temporal Pattern Analysis with an ADHD Case Study.*

## Overview

We adapt two algorithmic paradigms from stringology to discover recurrent temporal motifs in numerical EEG time series:

- **OPM** identifies subsequences that share the same rank ordering, capturing recurring ordinal trends (rises, falls, local extrema) independent of amplitude.
- **CTM** identifies subsequences whose min-Cartesian tree structure is identical, capturing hierarchical waveform morphology with greater tolerance to minor local variations.

Both methods are inherently amplitude-invariant and require no discretization or normalization prior to motif matching.

---

## Requirements

Install dependencies with:

```bash
pip install -r requirements.txt
```

Dependencies:
- Python 3.8+
- pandas
- numpy

---

## Data

This study uses the publicly available EEG dataset from:

> Nasrabadi, A.M., Allahverdy, A., Samavati, M., & Mohammadi, M. (2022). *EEG data for ADHD/control children.* IEEE DataPort.
> https://dx.doi.org/10.21227/rzfh-zn36

The dataset contains 19-channel EEG recordings from 61 children with ADHD and 60 typically developing controls (ages 7–12).

### Preparing the data

1. Download the dataset from the link above
2. Apply preprocessing (high-pass filter 0.5 Hz, low-pass filter 40 Hz, ICA artifact removal, 1-in-4 downsampling)
3. Save each participant's preprocessed EEG as a CSV file with 19 columns (one per channel)
4. Organize files into two folders:

```
data/
  control/    ← one CSV per control participant
  adhd/       ← one CSV per ADHD participant
```

A small synthetic sample dataset is provided in `data/sample/` for testing purposes.

---

## Usage

### OPM motif discovery

```bash
python opm_motif_discovery.py \
  --control data/control \
  --adhd data/adhd \
  --output results/opm_results.csv \
  --min_length 5 \
  --max_length 10 \
  --support 0.9
```

### CTM motif discovery

```bash
python ctm_motif_discovery.py \
  --control data/control \
  --adhd data/adhd \
  --output results/ctm_results.csv \
  --min_length 5 \
  --max_length 10 \
  --support 0.9
```

### Output format

Both scripts produce a CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `column` | EEG channel number (1–19) |
| `group` | Group label (`ADHD` or `Control`) |
| `motif` | Motif encoding (rank tuple for OPM, parentheses string for CTM) |
| `length` | Motif length in samples |
| `frequency` | Total number of occurrences across all sequences |
| `support` | Proportion of sequences containing the motif |
| `mean_value` | Position-wise mean amplitude across all motif instances |

---

## Repository structure

```
.
├── opm_motif_discovery.py    # OPM-based motif discovery
├── ctm_motif_discovery.py    # CTM-based motif discovery
├── requirements.txt          # Python dependencies
├── data/
│   └── sample/               # Synthetic sample data for testing
└── README.md
```

---

## Citation

If you use this code in your research, please cite:

```bibtex
@article{dahan2026stringology,
  author    = {Dahan, Anat and Ghazawi, Samah},
  title     = {Stringology-Based Motif Discovery for Electrophysiological
               Time Series: A Framework for Temporal Pattern Analysis
               with an ADHD Case Study},
  year      = {2026}
}
```

---

## License

This code is released under the MIT License.
