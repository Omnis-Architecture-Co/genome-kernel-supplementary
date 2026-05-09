import requests
import time
import os
import csv

ENSEMBL_SERVER = "https://rest.ensembl.org"
BASE_DIR = os.path.expanduser("~")
OUTPUT_DIR = os.path.join(BASE_DIR, "gene_sequences")
COMBINED_FILE = os.path.join(BASE_DIR, "all_100_genes_combined.txt")
MANIFEST_FILE = os.path.join(BASE_DIR, "gene_manifest.csv")
UPSTREAM_BP = 5000
SEPARATOR = ">>>>>><<<<<<<"

GENE_PANEL = {
    'Transcription': ['MYC', 'JUN', 'FOXP2', 'KLF4', 'IRF1'],
    'Chromatin': ['EZH2', 'HDAC1', 'KDM5B', 'SMARCA4', 'DNMT3A'],
    'Structural': ['COL1A1', 'LMNA', 'VIM', 'KRT18', 'FLNA'],
    'Cytoskeleton': ['ACTB', 'TUBB', 'MYH9', 'DYNC1H1', 'KIF5B'],
    'RNA_processing': ['SRSF1', 'HNRNPA1', 'DDX5', 'PRPF8', 'SF3B1'],
    'Cell_adhesion': ['CDH1', 'ITGB1', 'CD44', 'VCAM1', 'PECAM1'],
    'Cell_cycle': ['CDK2', 'CCND1', 'RB1', 'CDC20', 'PLK1'],
    'DNA_repair': ['BRCA1', 'ATM', 'MLH1', 'XRCC1', 'RAD51'],
    'Kinase': ['EGFR', 'BRAF', 'ABL1', 'JAK2', 'SRC'],
    'Ion_channel': ['SCN5A', 'KCNQ1', 'CFTR', 'CACNA1C', 'TRPV1'],
    'Phosphatase': ['PTEN', 'PTPN11', 'PPP2CA', 'PTPRC', 'DUSP1'],
    'Transport': ['SLC2A1', 'ABCB1', 'SLC6A3', 'ATP1A1', 'SLC12A2'],
    'Translation': ['EIF4E', 'EEF2', 'RPS6', 'EIF2AK2', 'MTOR'],
    'Apoptosis': ['BCL2', 'BAX', 'CASP3', 'XIAP', 'BID'],
    'Signaling': ['AKT1', 'KRAS', 'PIK3CA', 'NOTCH1', 'WNT3A'],
    'Immune': ['TNF', 'IL6', 'IFNG', 'CD4', 'TLR4'],
    'Proteolysis': ['MMP9', 'ADAM17', 'USP7', 'CTSD', 'CASP8'],
    'GTPase': ['RHOA', 'RAC1', 'CDC42', 'RAB7A', 'ARF1'],
    'Ubiquitin': ['UBE2I', 'RNF2', 'MDM2', 'NEDD4', 'TRIM28'],
    'Methylation': ['METTL3', 'PRMT1', 'SETD2', 'NSD1', 'DOT1L'],
}


def ensembl_lookup(symbol):
    ext = "/lookup/symbol/homo_sapiens/" + symbol
    r = requests.get(ENSEMBL_SERVER + ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        print("  WARNING: lookup failed for " + symbol)
        return None
    return r.json()


def fetch_sequence(chrom, start, end, strand):
    region = str(chrom) + ":" + str(start) + ".." + str(end) + ":" + str(strand)
    ext = "/sequence/region/human/" + region
    r = requests.get(ENSEMBL_SERVER + ext, headers={"Content-Type": "text/plain"})
    if not r.ok:
        print("  WARNING: sequence fetch failed for " + region)
        return None
    return r.text.strip()


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    manifest_rows = []
    combined_parts = []
    total = sum(len(g) for g in GENE_PANEL.values())
    count = 0
    failed = []
    char_offset = 0

    print("Fetching " + str(total) + " genes from Ensembl...")
    print("Output: " + OUTPUT_DIR)
    print("")

    for dept, genes in GENE_PANEL.items():
        for symbol in genes:
            count += 1
            print("[" + str(count) + "/" + str(total) + "] " + symbol + " (" + dept + ")... ", end="", flush=True)

            info = ensembl_lookup(symbol)
            time.sleep(0.15)

            if info is None:
                failed.append(symbol)
                print("FAILED (lookup)")
                continue

            chrom = info.get('seq_region_name', '')
            start = info.get('start')
            end = info.get('end')
            strand = info.get('strand')

            if strand == 1:
                full_start = max(1, start - UPSTREAM_BP)
                full_end = end
            else:
                full_start = start
                full_end = end + UPSTREAM_BP

            strand_char = "+" if strand == 1 else "-"
            print("chr" + str(chrom) + ":" + str(full_start) + "-" + str(full_end) + " (" + strand_char + ")... ", end="", flush=True)

            seq = fetch_sequence(chrom, full_start, full_end, strand)
            time.sleep(0.15)

            if seq is None:
                failed.append(symbol)
                print("FAILED (seq)")
                continue

            outfile = os.path.join(OUTPUT_DIR, symbol + "_" + dept + ".txt")
            with open(outfile, 'w') as f:
                f.write(seq)

            manifest_rows.append({
                'gene': symbol,
                'department': dept,
                'char_offset': char_offset,
                'seq_length': len(seq),
                'upstream_bp': UPSTREAM_BP,
                'upstream_byte_cutoff': UPSTREAM_BP // 4,
                'chromosome': "chr" + str(chrom),
                'strand': strand_char,
                'gene_start': start,
                'gene_end': end,
            })

            combined_parts.append(seq)
            char_offset += len(seq) + len(SEPARATOR)

            print("OK (" + str(len(seq)) + " bp)")

    combined_text = SEPARATOR.join(combined_parts)
    with open(COMBINED_FILE, 'w') as f:
        f.write(combined_text)
    print("")
    print("Combined file: " + COMBINED_FILE + " (" + str(len(combined_text)) + " chars)")

    if manifest_rows:
        with open(MANIFEST_FILE, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=manifest_rows[0].keys())
            writer.writeheader()
            writer.writerows(manifest_rows)
        print("Manifest: " + MANIFEST_FILE)

    print("")
    print("DONE: " + str(len(combined_parts)) + "/" + str(total) + " genes OK")
    if failed:
        print("FAILED (" + str(len(failed)) + "): " + ", ".join(failed))
    print("")
    print("Next: Feed all_100_genes_combined.txt through OMNIS")


if __name__ == "__main__":
    main()
