import os
import pandas as pd
from Bio import SeqIO
from datetime import datetime
import subprocess

# --------------------------
# Logger class
# --------------------------
class Logger:
    def __init__(self, log_file_path):
        self.log_file = open(log_file_path, 'w', encoding='utf-8')

    def write(self, message):
        print(message, end='')
        self.log_file.write(message)
        self.log_file.flush()

    def close(self):
        self.log_file.close()

# --------------------------
# Load top clusters
# --------------------------
def load_top_clusters(barcode_path, kmer):
    cluster_file = None
    for file in os.listdir(barcode_path):
        if file.endswith(f"_k{kmer}_clusters.csv") and file.lower().startswith(os.path.basename(barcode_path).lower()):
            cluster_file = os.path.join(barcode_path, file)
            break
    if not cluster_file:
        return None, None
    df = pd.read_csv(cluster_file)
    clusters = sorted(df['cluster'].unique())
    return df, clusters

# --------------------------
# Load sequences from FASTA
# --------------------------
def load_sequences(fasta_path, sequence_ids):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    matched_seqs = []
    for sid in sequence_ids:
        for rec_id, rec in seq_dict.items():
            if rec_id.startswith(sid):
                matched_seqs.append(rec)
                break
    return matched_seqs

# --------------------------
# Align sequences with MAFFT
# --------------------------
def align_with_mafft(input_fasta, output_fasta, mafft_path, logger):
    try:
        cmd = [mafft_path, "--auto", input_fasta]
        with open(output_fasta, "w") as out_handle:
            subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, check=True)
        logger.write(f"Aligned sequences saved to {output_fasta}\n")
    except subprocess.CalledProcessError as e:
        logger.write(f"Error aligning {input_fasta} with MAFFT:\n{e.stderr.decode()}\n")
    except FileNotFoundError:
        logger.write(f"MAFFT executable not found at: {mafft_path}\n")

# --------------------------
# Extract cluster sequences and save FASTA + align
# --------------------------
def extract_and_align_clusters(barcode_path, kmer, mafft_path, logger):
    output_dir = os.path.join(barcode_path, f"AlignedClusters-k{kmer}")
    os.makedirs(output_dir, exist_ok=True)

    logger.write(f"\n{'='*60}\nAnalyzing barcode: {barcode_path}, k={kmer}\n{'='*60}\n")

    df, clusters = load_top_clusters(barcode_path, kmer)
    if df is None:
        logger.write(f"No cluster CSV found for {os.path.basename(barcode_path)} with k={kmer}. Exiting.\n")
        return

    logger.write(f"Found clusters: {clusters}\n")

    fasta_file = os.path.join(barcode_path, f"{os.path.basename(barcode_path)}_combined.fasta")
    if not os.path.exists(fasta_file):
        logger.write(f"FASTA file {fasta_file} not found. Exiting.\n")
        return

    for cluster_id in clusters:
        cluster_fasta = os.path.join(output_dir, f"{os.path.basename(barcode_path)}_cluster{cluster_id}.fasta")
        aligned_fasta = os.path.join(output_dir, f"{os.path.basename(barcode_path)}_cluster{cluster_id}_aligned.fasta")

        # Skip if aligned file already exists
        if os.path.exists(aligned_fasta):
            logger.write(f"Skipping cluster {cluster_id}: alignment already exists at {aligned_fasta}\n")
            continue

        # If cluster fasta already exists but no alignment, reuse it
        if not os.path.exists(cluster_fasta):
            seq_ids = df[df['cluster'] == cluster_id]['sequence_id'].tolist()
            sequences = load_sequences(fasta_file, seq_ids)
            SeqIO.write(sequences, cluster_fasta, "fasta")
            logger.write(f"Saved {len(sequences)} sequences for cluster {cluster_id} to {cluster_fasta}\n")
        else:
            logger.write(f"Cluster fasta already exists for cluster {cluster_id}, using existing file {cluster_fasta}\n")

        # Align with MAFFT
        logger.write(f"Aligning cluster {cluster_id} with MAFFT...\n")
        align_with_mafft(cluster_fasta, aligned_fasta, mafft_path, logger)

# --------------------------
# Main function
# --------------------------
def main():
    base_path = r"C:\\Users\\ketgl\\OneDrive\\Desktop\\Sandbox\\Cactus2"
    mafft_path = r"C:\\Users\\ketgl\\OneDrive\\Documents\\GitHub\\Python-Programs\\BioTech\\Bioinformatics\\Tools\\MAFFT\\mafft-7.526-win64-signed (1)\\mafft-win\\mafft.bat"

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # List available barcodes
    all_barcodes = [f for f in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, f)) and f.lower().startswith('barcode')]
    print(f"Available barcodes:\n{', '.join(all_barcodes)}\n")

    barcode_choice = input("Enter barcode to analyze (number or full name): ").strip()
    if barcode_choice.isdigit():
        barcode_choice = f"barcode{barcode_choice}"

    if barcode_choice not in all_barcodes:
        print(f"Barcode {barcode_choice} not found. Exiting.")
        return

    kmer_input = input("Enter k-mer number to analyze (e.g., 6): ").strip()
    if not kmer_input.isdigit():
        print("Invalid k-mer input. Exiting.")
        return
    kmer = int(kmer_input)

    barcode_path = os.path.join(base_path, barcode_choice)

    log_file_path = os.path.join(barcode_path, f"barcode_clusters_alignment_k{kmer}_{timestamp}.log")
    logger = Logger(log_file_path)

    try:
        extract_and_align_clusters(barcode_path, kmer, mafft_path, logger)
    except Exception as e:
        logger.write(f"Unexpected error: {e}\n")
    finally:
        logger.close()
        print(f"\nLog file saved to: {log_file_path}")

if __name__ == "__main__":
    main()
