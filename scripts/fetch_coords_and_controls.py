"""
Fetch genomic coordinates for all 100 genes AND generate 100 random control windows.
Run: python fetch_coords_and_controls.py
Output: 
  - gene_panel_coordinates.csv (coordinates for Methods section)
  - random_control_sequences/ (100 random genomic windows)
  - all_random_controls_combined.txt (single file for OMNIS pipeline)
  - random_control_manifest.csv (position map)
"""

import requests
import time
import os
import csv
import random

ENSEMBL_SERVER = "https://rest.ensembl.org"
BASE_DIR = os.path.expanduser("~")
COORDS_FILE = os.path.join(BASE_DIR, "gene_panel_coordinates.csv")
RANDOM_DIR = os.path.join(BASE_DIR, "random_control_sequences")
RANDOM_COMBINED = os.path.join(BASE_DIR, "all_random_controls_combined.txt")
RANDOM_MANIFEST = os.path.join(BASE_DIR, "random_control_manifest.csv")
UPSTREAM_BP = 5000
SEPARATOR = "N" * 100

GENES = [
    'ABCB1', 'ABL1', 'ACTB', 'ADAM17', 'AKT1', 'ARF1', 'ATM', 'ATP1A1',
    'BAX', 'BCL2', 'BID', 'BRAF', 'BRCA1', 'CACNA1C', 'CASP3', 'CASP8',
    'CCND1', 'CD44', 'CD4', 'CDC20', 'CDC42', 'CDH1', 'CDK2', 'CFTR',
    'COL1A1', 'CTSD', 'DDX5', 'DNMT3A', 'DOT1L', 'DUSP1', 'DYNC1H1',
    'EEF2', 'EGFR', 'EIF2AK2', 'EIF4E', 'EZH2', 'FLNA', 'FOXP2', 'HDAC1',
    'HNRNPA1', 'IFNG', 'IL6', 'IRF1', 'ITGB1', 'JAK2', 'JUN', 'KCNQ1',
    'KDM5B', 'KIF5B', 'KLF4', 'KRAS', 'KRT18', 'LMNA', 'MDM2', 'METTL3',
    'MLH1', 'MMP9', 'MTOR', 'MYC', 'MYH9', 'NEDD4', 'NOTCH1', 'NSD1',
    'PECAM1', 'PIK3CA', 'PLK1', 'PPP2CA', 'PRMT1', 'PRPF8', 'PTEN',
    'PTPN11', 'PTPRC', 'RAB7A', 'RAC1', 'RAD51', 'RB1', 'RHOA', 'RNF2',
    'RPS6', 'SCN5A', 'SETD2', 'SF3B1', 'SLC12A2', 'SLC2A1', 'SLC6A3',
    'SMARCA4', 'SRC', 'SRSF1', 'TLR4', 'TNF', 'TRIM28', 'TRPV1', 'TUBB',
    'UBE2I', 'USP7', 'VCAM1', 'VIM', 'WNT3A', 'XIAP', 'XRCC1'
]

# Chromosome sizes for hg38 (approximate, for random window generation)
CHROM_SIZES = {
    '1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555,
    '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
    '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
    '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
    '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
    '21': 46709983, '22': 50818468, 'X': 156040895,
}


def ensembl_lookup(symbol):
    ext = "/lookup/symbol/homo_sapiens/" + symbol
    r = requests.get(ENSEMBL_SERVER + ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        return None
    return r.json()


def fetch_sequence(chrom, start, end, strand):
    region = str(chrom) + ":" + str(start) + ".." + str(end) + ":" + str(strand)
    ext = "/sequence/region/human/" + region
    r = requests.get(ENSEMBL_SERVER + ext, headers={"Content-Type": "text/plain"})
    if not r.ok:
        return None
    return r.text.strip()


def main():
    os.makedirs(RANDOM_DIR, exist_ok=True)

    # ============================================================
    # PART 1: Fetch coordinates for all 100 genes
    # ============================================================
    print("PART 1: Fetching gene coordinates...")
    coords = []
    gene_lengths = []

    for i, gene in enumerate(GENES):
        print("[" + str(i+1) + "/100] " + gene + "... ", end="", flush=True)
        info = ensembl_lookup(gene)
        time.sleep(0.15)

        if info is None:
            print("FAILED")
            continue

        chrom = info.get('seq_region_name', '')
        start = info.get('start')
        end = info.get('end')
        strand = info.get('strand')
        gene_length = end - start + 1

        coords.append({
            'gene': gene,
            'chromosome': "chr" + str(chrom),
            'gene_start': start,
            'gene_end': end,
            'strand': '+' if strand == 1 else '-',
            'gene_length_bp': gene_length,
            'upstream_bp': UPSTREAM_BP,
            'total_region_bp': gene_length + UPSTREAM_BP,
            'ensembl_id': info.get('id', ''),
        })
        gene_lengths.append(gene_length + UPSTREAM_BP)
        print("OK (" + str(gene_length) + " bp)")

    with open(COORDS_FILE, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=coords[0].keys())
        writer.writeheader()
        writer.writerows(coords)
    print("Coordinates saved: " + COORDS_FILE)

    # ============================================================
    # PART 2: Generate 100 random genomic windows as controls
    # ============================================================
    print("\nPART 2: Generating 100 random control windows...")

    # Match the length distribution of real genes
    # For each real gene, pick a random intergenic window of same length
    random.seed(42)
    chroms = list(CHROM_SIZES.keys())

    random_parts = []
    random_manifest = []
    char_offset = 0

    for i in range(100):
        # Pick a random chromosome and position
        chrom = random.choice(chroms)
        chrom_size = CHROM_SIZES[chrom]

        # Match length to corresponding real gene
        target_length = gene_lengths[i] if i < len(gene_lengths) else 20000
        
        # Pick random start, avoiding edges
        max_start = chrom_size - target_length - 10000
        if max_start < 10000:
            max_start = 10000
        rand_start = random.randint(10000, max_start)
        rand_end = rand_start + target_length

        # Random strand
        strand = random.choice([1, -1])

        print("[" + str(i+1) + "/100] chr" + chrom + ":" + str(rand_start) + "-" + str(rand_end) + "... ", end="", flush=True)

        seq = fetch_sequence(chrom, rand_start, rand_end, strand)
        time.sleep(0.15)

        if seq is None:
            print("FAILED")
            continue

        outfile = os.path.join(RANDOM_DIR, "random_" + str(i+1).zfill(3) + "_chr" + chrom + ".txt")
        with open(outfile, 'w') as f:
            f.write(seq)

        random_manifest.append({
            'window_id': "random_" + str(i+1).zfill(3),
            'chromosome': "chr" + chrom,
            'start': rand_start,
            'end': rand_end,
            'strand': '+' if strand == 1 else '-',
            'length_bp': len(seq),
            'char_offset': char_offset,
            'seq_length': len(seq),
            'upstream_byte_cutoff': UPSTREAM_BP // 4,
        })

        random_parts.append(seq)
        char_offset += len(seq) + len(SEPARATOR)
        print("OK (" + str(len(seq)) + " bp)")

    # Write combined random file
    combined = SEPARATOR.join(random_parts)
    with open(RANDOM_COMBINED, 'w') as f:
        f.write(combined)
    print("\nCombined random file: " + RANDOM_COMBINED + " (" + str(len(combined)) + " chars)")

    # Write random manifest
    if random_manifest:
        with open(RANDOM_MANIFEST, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=random_manifest[0].keys())
            writer.writeheader()
            writer.writerows(random_manifest)
        print("Random manifest: " + RANDOM_MANIFEST)

    print("\n" + "="*60)
    print("DONE")
    print("1. Gene coordinates: " + COORDS_FILE)
    print("2. Random controls: " + RANDOM_COMBINED)
    print("3. Random manifest: " + RANDOM_MANIFEST)
    print("\nNext: Run all_random_controls_combined.txt through OMNIS")
    print("Then upload the token CSV + random_control_manifest.csv")


if __name__ == "__main__":
    main()
