"""
OMNIS Permutation Test for Upstream PIC Enrichment
====================================================
Tests whether upstream tokens carry PIC-associated carrier proteins
at a rate exceeding chance expectation.

Randomly samples tokens from the combined pool and counts PIC tokens.
Repeats 10,000 times per seed. Run with multiple seeds to confirm stability.

Inputs:
  - omnis_100gene_all_tokens.csv (full token export)
  - gene_manifest.csv (byte offsets)
  - vocabulary_human_1932words.csv (vocabulary with carrier genes)

Outputs:
  - permutation_results_seed_XX.txt (results for each seed)
  - permutation_null_distribution_seed_XX.csv (null values for plotting)

Usage:
  python permutation_test.py              (runs all 3 seeds)
  python permutation_test.py 42           (runs single seed)
  python permutation_test.py 42 123 456   (runs specified seeds)
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import sys

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TOKENS_FILE = os.path.join(BASE_DIR, "omnis_100gene_all_tokens.csv")
MANIFEST_FILE = os.path.join(BASE_DIR, "gene_manifest.csv")
VOCAB_FILE = os.path.join(BASE_DIR, "vocabulary_human_1932words.csv")
N_PERMUTATIONS = 10000

PIC_UPPER = {
    'TBP', 'TAF1', 'TAF2', 'TAF3', 'TAF4', 'TAF5', 'TAF6', 'TAF7',
    'TAF8', 'TAF9', 'TAF10', 'TAF11', 'TAF12', 'TAF13',
    'GTF2A1', 'GTF2A2', 'GTF2B',
    'POLR2A', 'POLR2B', 'POLR2C', 'POLR2D', 'POLR2E', 'POLR2F',
    'POLR2G', 'POLR2H', 'POLR2I', 'POLR2J', 'POLR2K', 'POLR2L',
    'GTF2F1', 'GTF2F2', 'GTF2E1', 'GTF2E2',
    'GTF2H1', 'GTF2H2', 'GTF2H3', 'GTF2H4', 'GTF2H5',
    'ERCC2', 'ERCC3', 'CDK7', 'CCNH', 'MNAT1',
    'MED1', 'MED4', 'MED6', 'MED7', 'MED8', 'MED9', 'MED10', 'MED11',
    'MED12', 'MED13', 'MED14', 'MED15', 'MED16', 'MED17', 'MED18',
    'MED19', 'MED20', 'MED21', 'MED22', 'MED23', 'MED24', 'MED25',
    'MED26', 'MED27', 'MED28', 'MED29', 'MED30', 'MED31',
}


def run_permutation(seed):
    print("=" * 60)
    print("PERMUTATION TEST (seed=" + str(seed) + ", n=" + str(N_PERMUTATIONS) + ")")
    print("=" * 60)

    # Load data
    tokens = pd.read_csv(TOKENS_FILE)
    manifest = pd.read_csv(MANIFEST_FILE)
    vocab = pd.read_csv(VOCAB_FILE)

    tokens['hex_norm'] = tokens['hex'].str.replace('"', '').str.strip().str.upper()
    vocab['hex_norm'] = vocab['word_hex'].str.strip().str.upper()

    manifest['byte_start'] = manifest['char_offset'] // 4
    manifest['byte_end'] = (manifest['char_offset'] + manifest['seq_length']) // 4
    manifest['upstream_end'] = manifest['byte_start'] + manifest['upstream_byte_cutoff']

    # Build PIC token set
    pic_tokens = set()
    for _, row in vocab.iterrows():
        carriers = str(row.get('all_carrier_genes', '')).split(';')
        carriers = [c.strip() for c in carriers if c.strip()]
        if any(c.upper() in PIC_UPPER for c in carriers):
            pic_tokens.add(row['hex_norm'])

    # Assign tokens to upstream vs genebody
    manifest_sorted = manifest.sort_values('byte_start').reset_index(drop=True)
    m_starts = manifest_sorted['byte_start'].values
    m_ends = manifest_sorted['byte_end'].values

    is_upstream = []
    is_assigned = []
    for pos in tokens['position']:
        idx = np.searchsorted(m_starts, pos, side='right') - 1
        if 0 <= idx < len(manifest_sorted) and pos < m_ends[idx]:
            row = manifest_sorted.iloc[idx]
            is_assigned.append(True)
            is_upstream.append(pos < row['upstream_end'])
        else:
            is_assigned.append(False)
            is_upstream.append(False)

    tokens['assigned'] = is_assigned
    tokens['is_upstream'] = is_upstream

    assigned = tokens[tokens['assigned']]
    upstream_hexes = assigned[assigned['is_upstream']]['hex_norm'].values
    all_hexes = assigned['hex_norm'].values

    # Observed
    observed_pic = sum(1 for h in upstream_hexes if h in pic_tokens)
    n_upstream = len(upstream_hexes)

    print("Upstream tokens: " + str(n_upstream))
    print("Observed PIC tokens: " + str(observed_pic))
    print("Observed PIC rate: " + str(round(observed_pic / n_upstream * 100, 2)) + "%")
    print("")
    print("Running " + str(N_PERMUTATIONS) + " permutations...")

    # Permutation
    np.random.seed(seed)
    null_pics = np.zeros(N_PERMUTATIONS, dtype=int)

    for i in range(N_PERMUTATIONS):
        sample = np.random.choice(all_hexes, size=n_upstream, replace=False)
        null_pics[i] = sum(1 for h in sample if h in pic_tokens)

        if (i + 1) % 2000 == 0:
            print("  " + str(i + 1) + " / " + str(N_PERMUTATIONS) + " complete")

    null_mean = null_pics.mean()
    null_std = null_pics.std()
    z_score = (observed_pic - null_mean) / null_std
    p_value = (np.sum(null_pics >= observed_pic) + 1) / (N_PERMUTATIONS + 1)

    # Results
    lines = []
    lines.append("=" * 60)
    lines.append("PERMUTATION TEST RESULTS")
    lines.append("=" * 60)
    lines.append("Seed: " + str(seed))
    lines.append("Permutations: " + str(N_PERMUTATIONS))
    lines.append("Upstream tokens: " + str(n_upstream))
    lines.append("Total pool tokens: " + str(len(all_hexes)))
    lines.append("")
    lines.append("Observed PIC tokens: " + str(observed_pic))
    lines.append("Observed PIC rate: " + str(round(observed_pic / n_upstream * 100, 2)) + "%")
    lines.append("")
    lines.append("Null mean: " + str(round(null_mean, 1)))
    lines.append("Null std: " + str(round(null_std, 1)))
    lines.append("Null PIC rate: " + str(round(null_mean / n_upstream * 100, 2)) + "%")
    lines.append("")
    lines.append("z-score: " + str(round(z_score, 3)))
    lines.append("p-value: " + str(round(p_value, 6)))
    lines.append("")
    lines.append("Null distribution range: " + str(null_pics.min()) + " - " + str(null_pics.max()))
    lines.append("Null 95th percentile: " + str(round(np.percentile(null_pics, 95), 1)))
    lines.append("Null 99th percentile: " + str(round(np.percentile(null_pics, 99), 1)))

    result_text = '\n'.join(lines)
    print("")
    print(result_text)

    # Save results
    result_file = os.path.join(BASE_DIR, "permutation_results_seed_" + str(seed) + ".txt")
    with open(result_file, 'w') as f:
        f.write(result_text)
    print("Results saved: " + result_file)

    # Save null distribution for plotting
    null_file = os.path.join(BASE_DIR, "permutation_null_seed_" + str(seed) + ".csv")
    pd.DataFrame({'null_pic_count': null_pics}).to_csv(null_file, index=False)
    print("Null distribution saved: " + null_file)

    return z_score, p_value


def main():
    if len(sys.argv) > 1:
        seeds = [int(s) for s in sys.argv[1:]]
    else:
        seeds = [42, 123, 456]

    results = []
    for seed in seeds:
        z, p = run_permutation(seed)
        results.append((seed, z, p))
        print("")

    if len(results) > 1:
        print("=" * 60)
        print("CROSS-SEED SUMMARY")
        print("=" * 60)
        for seed, z, p in results:
            print("  Seed " + str(seed) + ": z = " + str(round(z, 3)) + ", p = " + str(round(p, 6)))

        z_values = [r[1] for r in results]
        p_values = [r[2] for r in results]
        print("")
        print("  Mean z: " + str(round(np.mean(z_values), 3)))
        print("  z range: " + str(round(min(z_values), 3)) + " - " + str(round(max(z_values), 3)))
        print("  All p < 0.05: " + str(all(p < 0.05 for p in p_values)))

        summary_file = os.path.join(BASE_DIR, "permutation_cross_seed_summary.txt")
        with open(summary_file, 'w') as f:
            f.write("Cross-seed permutation test summary\n")
            f.write("Seeds: " + str(seeds) + "\n")
            for seed, z, p in results:
                f.write("Seed " + str(seed) + ": z = " + str(round(z, 3)) + ", p = " + str(round(p, 6)) + "\n")
            f.write("Mean z: " + str(round(np.mean(z_values), 3)) + "\n")
            f.write("All significant: " + str(all(p < 0.05 for p in p_values)) + "\n")
        print("Summary saved: " + summary_file)


if __name__ == "__main__":
    main()
