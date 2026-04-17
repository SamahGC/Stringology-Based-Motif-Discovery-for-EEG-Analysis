"""
CTM-based motif discovery for EEG time series.

Cartesian Tree Matching (CTM) identifies recurrent subsequences whose
min-Cartesian tree structure is identical, capturing hierarchical waveform
morphology independent of absolute amplitude values.

Usage:
    python ctm_motif_discovery.py \
        --control data/control \
        --adhd data/adhd \
        --output results/ctm_results.csv \
        --min_length 5 \
        --max_length 10 \
        --support 0.9
"""

import os
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict


# =========================
#   CARTESIAN TREE ENCODING
# =========================
def cartesian_tree_encoding_numeric(subseq):
    """Return compact parentheses shape of min-Cartesian tree for numeric subsequence."""
    stack = []
    left = {}
    right = {}

    for i, val in enumerate(subseq):
        last = None
        while stack and subseq[stack[-1]] > val:
            last = stack.pop()
        if stack:
            right[stack[-1]] = i
        if last is not None:
            left[i] = last
        stack.append(i)

    root = stack[0] if stack else None

    def encode(node):
        if node is None:
            return ""
        return "(" + encode(left.get(node)) + encode(right.get(node)) + ")"

    return encode(root)


# =========================
#   MOTIF DISCOVERY CLASS
# =========================
class MotifDiscoveryCTM:
    def __init__(self, min_length=5, max_length=10):
        self.min_length = min_length
        self.max_length = max_length

    def extract_sequences_from_csv_per_column(self, csv_path, expected_columns=19):
        """Extract sequences per column from one CSV file."""
        try:
            df = pd.read_csv(csv_path)
            column_sequences = {}
            for idx in range(expected_columns):
                col_number = idx + 1
                if idx < len(df.columns):
                    seq = df.iloc[:, idx].dropna().astype(float).tolist()
                else:
                    seq = []
                column_sequences[col_number] = seq
            return column_sequences
        except Exception as e:
            print(f"Warning: error reading {csv_path}: {e}")
            return {}

    def find_motifs_single_length(self, sequences, motif_length, min_support=0.9):
        """Find Cartesian Tree motifs of fixed length with support filtering."""
        motif_counts = defaultdict(int)
        motif_sequences = defaultdict(set)
        motif_real_values = defaultdict(list)

        for seq_idx, numeric_seq in enumerate(sequences):
            n = len(numeric_seq)
            if n < motif_length:
                continue
            for i in range(n - motif_length + 1):
                substring = numeric_seq[i:i + motif_length]
                motif = cartesian_tree_encoding_numeric(substring)
                motif_counts[motif] += 1
                motif_sequences[motif].add(seq_idx)
                motif_real_values[motif].append(substring)

        motifs = {}
        for motif, count in motif_counts.items():
            support = len(motif_sequences[motif]) / len(sequences)
            if support >= min_support:
                subseqs = np.array(motif_real_values[motif])
                mean_sequence = np.mean(subseqs, axis=0).round(4).tolist()
                motifs[motif] = {
                    'frequency': count,
                    'support': round(support, 3),
                    'mean_value': mean_sequence
                }
        return motifs


# =========================
#   REMOVE SUBSTRING MOTIFS
# =========================
def remove_substring_motifs(motifs_dict):
    """Remove motifs fully contained in another motif (keep maximal motifs)."""
    motifs_sorted = sorted(motifs_dict.keys(), key=lambda x: len(x), reverse=True)
    maximal_motifs = {}
    for m in motifs_sorted:
        is_substring = False
        for kept in maximal_motifs:
            if m in kept:
                is_substring = True
                break
        if not is_substring:
            maximal_motifs[m] = motifs_dict[m]
    return maximal_motifs


# =========================
#   DATA LOADING
# =========================
def process_group_per_column(input_folder, motif_finder, expected_columns=19):
    """Load sequences per column from all CSV files in a folder."""
    files = sorted([f for f in os.listdir(input_folder) if f.endswith(".csv")])
    all_column_sequences = defaultdict(list)
    for file_name in files:
        col_seqs = motif_finder.extract_sequences_from_csv_per_column(
            os.path.join(input_folder, file_name),
            expected_columns=expected_columns
        )
        for col, seq in col_seqs.items():
            all_column_sequences[col].append(seq)
    print(f"Loaded {len(files)} files from {input_folder}")
    return all_column_sequences


# =========================
#   MOTIF COMPARISON AND SAVE
# =========================
def compare_and_save(control_columns, adhd_columns, motif_finder, output_file, min_support=0.9):
    rows = []
    all_columns = set(control_columns.keys()) | set(adhd_columns.keys())

    for col in sorted(all_columns):
        ctrl_seqs = control_columns.get(col, [])
        adhd_seqs = adhd_columns.get(col, [])

        for L in range(motif_finder.min_length, motif_finder.max_length + 1):
            ctrl_motifs = motif_finder.find_motifs_single_length(ctrl_seqs, L, min_support)
            adhd_motifs = motif_finder.find_motifs_single_length(adhd_seqs, L, min_support)

            ctrl_motifs = remove_substring_motifs(ctrl_motifs)
            adhd_motifs = remove_substring_motifs(adhd_motifs)

            ctrl_unique = set(ctrl_motifs.keys()) - set(adhd_motifs.keys())
            adhd_unique = set(adhd_motifs.keys()) - set(ctrl_motifs.keys())

            for motif in ctrl_unique:
                stats = ctrl_motifs[motif]
                rows.append({
                    'column': col, 'group': 'Control',
                    'motif': motif,
                    'length': L,
                    'frequency': stats['frequency'],
                    'support': stats['support'],
                    'mean_value': str(stats['mean_value'])
                })

            for motif in adhd_unique:
                stats = adhd_motifs[motif]
                rows.append({
                    'column': col, 'group': 'ADHD',
                    'motif': motif,
                    'length': L,
                    'frequency': stats['frequency'],
                    'support': stats['support'],
                    'mean_value': str(stats['mean_value'])
                })

        print(f"Finished channel {col}")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    pd.DataFrame(rows).to_csv(output_file, index=False)
    print(f"CTM results saved to {output_file}")


# =========================
#   MAIN
# =========================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CTM-based EEG motif discovery")
    parser.add_argument("--control", required=True, help="Path to control group CSV folder")
    parser.add_argument("--adhd", required=True, help="Path to ADHD group CSV folder")
    parser.add_argument("--output", required=True, help="Output CSV file path")
    parser.add_argument("--min_length", type=int, default=5, help="Minimum motif length")
    parser.add_argument("--max_length", type=int, default=10, help="Maximum motif length")
    parser.add_argument("--support", type=float, default=0.9, help="Minimum support threshold")
    args = parser.parse_args()

    motif_finder = MotifDiscoveryCTM(min_length=args.min_length, max_length=args.max_length)
    control_columns = process_group_per_column(args.control, motif_finder)
    adhd_columns = process_group_per_column(args.adhd, motif_finder)
    compare_and_save(control_columns, adhd_columns, motif_finder, args.output, args.support)
