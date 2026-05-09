"""
OMNIS Upstream PIC Analysis
============================
Maps vocabulary tokens to genes, identifies PIC carrier proteins
in upstream vs gene body regions, computes paired t-test and effect size.

Inputs:
  - omnis_100gene_all_tokens.csv (full token export from OMNIS pipeline)
  - gene_manifest.csv (byte offsets for each gene in combined file)
  - vocabulary_human_1932words.csv (vocabulary with carrier genes)

Outputs:
  - 100gene_PIC_analysis.csv (per-gene results)
  - upstream_pic_summary.txt (summary statistics)

Usage: python upstream_pic_analysis.py
"""

import pandas as pd
import numpy as np
from scipy import stats
from collections import Counter
import os
import sys

# Configuration
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TOKENS_FILE = os.path.join(BASE_DIR, "omnis_100gene_all_tokens.csv")
MANIFEST_FILE = os.path.join(BASE_DIR, "gene_manifest.csv")
VOCAB_FILE = os.path.join(BASE_DIR, "vocabulary_human_1932words.csv")
OUTPUT_CSV = os.path.join(BASE_DIR, "100gene_PIC_analysis.csv")
OUTPUT_SUMMARY = os.path.join(BASE_DIR, "upstream_pic_summary.txt")

# 75 PIC components across 7 subcomplexes
PIC_PROTEINS = {
    # TFIID/TBP (14)
    'TBP', 'TAF1', 'TAF2', 'TAF3', 'TAF4', 'TAF5', 'TAF6', 'TAF7',
    'TAF8', 'TAF9', 'TAF10', 'TAF11', 'TAF12', 'TAF13',
    # TFIIA (2)
    'GTF2A1', 'GTF2A2',
    # TFIIB (1)
    'GTF2B',
    # RNA Pol II + TFIIF (14)
    'POLR2A', 'POLR2B', 'POLR2C', 'POLR2D', 'POLR2E', 'POLR2F',
    'POLR2G', 'POLR2H', 'POLR2I', 'POLR2J', 'POLR2K', 'POLR2L',
    'GTF2F1', 'GTF2F2',
    # TFIIE (2)
    'GTF2E1', 'GTF2E2',
    # TFIIH (10)
    'GTF2H1', 'GTF2H2', 'GTF2H3', 'GTF2H4', 'GTF2H5',
    'ERCC2', 'ERCC3', 'CDK7', 'CCNH', 'MNAT1',
    # Mediator (31)
    'MED1', 'MED4', 'MED6', 'MED7', 'MED8', 'MED9', 'MED10', 'MED11',
    'MED12', 'MED13', 'MED14', 'MED15', 'MED16', 'MED17', 'MED18',
    'MED19', 'MED20', 'MED21', 'MED22', 'MED23', 'MED24', 'MED25',
    'MED26', 'MED27', 'MED28', 'MED29', 'MED30', 'MED31',
}
PIC_UPPER = {p.upper() for p in PIC_PROTEINS}


def load_data():
    print("Loading data...")
    tokens = pd.read_csv(TOKENS_FILE)
    manifest = pd.read_csv(MANIFEST_FILE)
    vocab = pd.read_csv(VOCAB_FILE)

    tokens['hex_norm'] = tokens['hex'].str.replace('"', '').str.strip().str.upper()
    vocab['hex_norm'] = vocab['word_hex'].str.strip().str.upper()

    manifest['byte_start'] = manifest['char_offset'] // 4
    manifest['byte_end'] = (manifest['char_offset'] + manifest['seq_length']) // 4
    manifest['upstream_end'] = manifest['byte_start'] + manifest['upstream_byte_cutoff']

    return tokens, manifest, vocab


def assign_tokens_to_genes(tokens, manifest):
    print("Assigning tokens to genes...")
    manifest_sorted = manifest.sort_values('byte_start').reset_index(drop=True)
    m_starts = manifest_sorted['byte_start'].values
    m_ends = manifest_sorted['byte_end'].values

    gene_list = []
    region_list = []

    for pos in tokens['position']:
        idx = np.searchsorted(m_starts, pos, side='right') - 1
        if 0 <= idx < len(manifest_sorted) and pos < m_ends[idx]:
            row = manifest_sorted.iloc[idx]
            gene_list.append(row['gene'])
            region_list.append('upstream' if pos < row['upstream_end'] else 'genebody')
        else:
            gene_list.append(None)
            region_list.append(None)

    tokens['gene'] = gene_list
    tokens['region'] = region_list
    return tokens, manifest_sorted


def build_pic_token_set(vocab):
    print("Building PIC token set...")
    pic_tokens = set()
    token_pic_proteins = {}

    for _, row in vocab.iterrows():
        carriers = str(row.get('all_carrier_genes', '')).split(';')
        carriers = [c.strip() for c in carriers if c.strip()]
        pic_hits = [c for c in carriers if c.upper() in PIC_UPPER]
        if pic_hits:
            pic_tokens.add(row['hex_norm'])
            token_pic_proteins[row['hex_norm']] = pic_hits

    print("  Vocabulary tokens carrying PIC: " + str(len(pic_tokens)) + " / " + str(len(vocab)))
    print("  Baseline rate: " + str(round(len(pic_tokens) / len(vocab) * 100, 1)) + "%")
    return pic_tokens, token_pic_proteins


def analyze_genes(tokens, manifest_sorted, pic_tokens, token_pic_proteins):
    print("Analyzing per-gene PIC enrichment...")
    assigned = tokens[tokens['gene'].notna()]
    upstream = assigned[assigned['region'] == 'upstream']
    genebody = assigned[assigned['region'] == 'genebody']

    results = []
    for _, mrow in manifest_sorted.iterrows():
        gene = mrow['gene']
        dept = mrow['department']

        gene_up = upstream[upstream['gene'] == gene]
        gene_gb = genebody[genebody['gene'] == gene]

        up_pic = sum(1 for h in gene_up['hex_norm'] if h in pic_tokens)
        gb_pic = sum(1 for h in gene_gb['hex_norm'] if h in pic_tokens)

        up_pct = (up_pic / len(gene_up) * 100) if len(gene_up) > 0 else 0
        gb_pct = (gb_pic / len(gene_gb) * 100) if len(gene_gb) > 0 else 0

        # Get specific PIC proteins found upstream
        up_proteins = set()
        for h in gene_up['hex_norm']:
            if h in token_pic_proteins:
                up_proteins.update(token_pic_proteins[h])

        results.append({
            'gene': gene,
            'department': dept,
            'upstream_tokens': len(gene_up),
            'upstream_pic_tokens': up_pic,
            'upstream_pic_pct': round(up_pct, 2),
            'genebody_tokens': len(gene_gb),
            'genebody_pic_tokens': gb_pic,
            'genebody_pic_pct': round(gb_pct, 2),
            'upstream_pic_proteins': ';'.join(sorted(up_proteins)),
        })

    return pd.DataFrame(results), upstream, genebody


def compute_statistics(df, upstream, genebody, pic_tokens, tokens):
    lines = []
    lines.append("=" * 70)
    lines.append("OMNIS UPSTREAM PIC ANALYSIS RESULTS")
    lines.append("=" * 70)

    assigned = tokens[tokens['gene'].notna()]
    lines.append("")
    lines.append("Token counts:")
    lines.append("  Total tokens: " + str(len(tokens)))
    lines.append("  Assigned to genes: " + str(len(assigned)))
    lines.append("  Upstream: " + str(len(upstream)))
    lines.append("  Gene body: " + str(len(genebody)))

    genes_with_pic = (df['upstream_pic_tokens'] > 0).sum()
    lines.append("")
    lines.append("PIC enrichment:")
    lines.append("  Genes with upstream PIC carriers: " + str(genes_with_pic) + " / " + str(len(df)))

    up_pcts = df['upstream_pic_pct'].values
    gb_pcts = df['genebody_pic_pct'].values

    lines.append("  Mean upstream PIC %: " + str(round(up_pcts.mean(), 2)))
    lines.append("  Mean genebody PIC %: " + str(round(gb_pcts.mean(), 2)))

    t, p = stats.ttest_rel(up_pcts, gb_pcts)
    diff = up_pcts - gb_pcts
    cohens_d = diff.mean() / diff.std()

    lines.append("")
    lines.append("Paired t-test (upstream vs genebody):")
    lines.append("  t = " + str(round(t, 3)))
    lines.append("  p = " + str(round(p, 6)))
    lines.append("  Cohen's d = " + str(round(cohens_d, 3)))

    # Vocabulary baseline
    baseline = len(pic_tokens) / 1932 * 100
    up_rate = up_pcts.mean()
    enrichment = up_rate / baseline

    lines.append("")
    lines.append("Baseline comparison:")
    lines.append("  Vocabulary baseline PIC rate: " + str(round(baseline, 1)) + "%")
    lines.append("  Upstream enrichment vs baseline: " + str(round(enrichment, 2)) + "x")

    # PIC protein prevalence
    all_upstream_pic = Counter()
    for _, row in df.iterrows():
        if row['upstream_pic_proteins']:
            for p_name in row['upstream_pic_proteins'].split(';'):
                if p_name:
                    all_upstream_pic[p_name] += 1

    lines.append("")
    lines.append("Top PIC proteins found upstream (gene count):")
    for protein, count in all_upstream_pic.most_common(20):
        lines.append("  " + protein + ": " + str(count) + " / " + str(len(df)) + " genes")

    # Four negative genes
    no_pic = df[df['upstream_pic_tokens'] == 0]
    lines.append("")
    lines.append("Genes with NO upstream PIC: " + str(len(no_pic)))
    for _, row in no_pic.iterrows():
        lines.append("  " + row['gene'] + " (" + row['department'] + "): " + str(row['upstream_tokens']) + " upstream tokens")

    return '\n'.join(lines)


def main():
    tokens, manifest, vocab = load_data()
    tokens, manifest_sorted = assign_tokens_to_genes(tokens, manifest)
    pic_tokens, token_pic_proteins = build_pic_token_set(vocab)
    df, upstream, genebody = analyze_genes(tokens, manifest_sorted, pic_tokens, token_pic_proteins)

    # Save per-gene results
    df.to_csv(OUTPUT_CSV, index=False)
    print("Per-gene results saved: " + OUTPUT_CSV)

    # Compute and save summary
    summary = compute_statistics(df, upstream, genebody, pic_tokens, tokens)
    with open(OUTPUT_SUMMARY, 'w') as f:
        f.write(summary)
    print("Summary saved: " + OUTPUT_SUMMARY)
    print("")
    print(summary)


if __name__ == "__main__":
    main()
