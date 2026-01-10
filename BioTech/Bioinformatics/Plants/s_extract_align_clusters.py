import os
import re
import pickle
import pandas as pd
from Bio import SeqIO
from datetime import datetime
import subprocess
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import time
import signal
import sys

# --------------------------
# Logger class
# --------------------------
class BufferedLogger:
    def __init__(self, log_file_path, buffer_size=1000):
        self.log_file = open(log_file_path, 'a', encoding='utf-8')  # append mode
        self.buffer = []
        self.buffer_size = buffer_size

    def write(self, message):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        full_msg = f"[{timestamp}] {message}\n"
        print(full_msg, end='')
        self.buffer.append(full_msg)
        if len(self.buffer) >= self.buffer_size:
            self.flush()
        else:
            self.flush()  # immediate flush to disk to protect against termination

    def flush(self):
        if self.buffer:
            self.log_file.write(''.join(self.buffer))
            self.log_file.flush()
            os.fsync(self.log_file.fileno())  # ensure written to disk
            self.buffer = []

    def close(self):
        self.flush()
        self.log_file.close()


# --------------------------
# Normalize FASTA record IDs
# --------------------------
def normalize_id(record_id: str) -> str:
    match = re.match(r'([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})', record_id)
    return match.group(1) if match else record_id

def clean_sequence_id(seq_id: str) -> str:
    match = re.match(r'([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})', seq_id)
    return match.group(1) if match else seq_id


# --------------------------
# Load top clusters
# --------------------------
def load_top_clusters(barcode_path, kmer):
    barcode_name = Path(barcode_path).name.lower()
    cluster_files = [
        f for f in os.listdir(barcode_path)
        if f.endswith(f"_k{kmer}_clusters.csv") and f.lower().startswith(barcode_name)
    ]

    if not cluster_files:
        return None, None

    cluster_file = os.path.join(barcode_path, cluster_files[0])
    try:
        df = pd.read_csv(cluster_file, usecols=['sequence_id', 'cluster'], dtype={'cluster': 'category'})
    except Exception:
        df = pd.read_csv(cluster_file)

    df['sequence_id'] = df['sequence_id'].astype(str).map(clean_sequence_id)
    clusters = sorted(df['cluster'].unique())
    return df, clusters


# --------------------------
# Load sequences
# --------------------------
def load_sequences_optimized(fasta_path, sequence_ids):
    sequence_ids_set = set(map(str, sequence_ids))
    matched_seqs = []
    with open(fasta_path, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            norm_id = normalize_id(record.id)
            if norm_id in sequence_ids_set:
                matched_seqs.append(record)
                sequence_ids_set.remove(norm_id)
            if not sequence_ids_set:
                break

    if sequence_ids_set:
        print(f"⚠️ Warning: {len(sequence_ids_set)} IDs not found in FASTA: {list(sequence_ids_set)[:10]}...")
    return matched_seqs


# --------------------------
# Create/load FASTA index
# --------------------------
def create_or_load_fasta_index(fasta_path, index_path, logger):
    if os.path.exists(index_path):
        logger.write(f"Loading existing FASTA index from {index_path}")
        with open(index_path, "rb") as f:
            fasta_index = pickle.load(f)
        logger.write(f"Loaded FASTA index with {len(fasta_index)} entries")
    else:
        logger.write(f"Creating new FASTA index at {index_path}")
        fasta_index = {}
        with open(fasta_path, "r") as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                norm_id = normalize_id(record.id)
                fasta_index[norm_id] = record
        with open(index_path, "wb") as f:
            pickle.dump(fasta_index, f)
        logger.write(f"FASTA index created with {len(fasta_index)} entries")
    return fasta_index


def load_sequences_with_index(fasta_index, sequence_ids):
    matched_seqs = []
    for sid in sequence_ids:
        if sid in fasta_index:
            matched_seqs.append(fasta_index[sid])
        else:
            print(f"⚠️ Warning: ID {sid} not found in FASTA index")
    return matched_seqs


# --------------------------
# MAFFT alignment
# --------------------------
def align_with_mafft_parallel(input_fasta, output_fasta, mafft_path, logger, threads=2):
    temp_output = output_fasta + ".tmp"
    try:
        cmd = [mafft_path, "--auto", "--thread", str(threads), input_fasta]
        with open(temp_output, "w") as out_handle:
            subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, check=True)
        os.replace(temp_output, output_fasta)
        logger.write(f"Aligned sequences saved to {output_fasta} using {threads} threads")
        return True
    except subprocess.CalledProcessError as e:
        logger.write(f"Error aligning {input_fasta} with MAFFT:\n{e.stderr.decode()}")
        if os.path.exists(temp_output):
            os.remove(temp_output)
        return False
    except FileNotFoundError:
        logger.write(f"MAFFT executable not found at: {mafft_path}")
        return False


# --------------------------
# Process single cluster
# --------------------------
def process_cluster(args):
    cluster_id, df_cluster, fasta_path, output_dir, barcode_name, mafft_path, use_index, fasta_index = args
    start_time = time.time()
    print(f"\nProcessing cluster {cluster_id}... Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    cluster_fasta = os.path.join(output_dir, f"{barcode_name}_cluster{cluster_id}.fasta")
    aligned_fasta = os.path.join(output_dir, f"{barcode_name}_cluster{cluster_id}_aligned.fasta")

    if os.path.exists(aligned_fasta):
        print(f"Skipping cluster {cluster_id}: alignment already exists")
        return f"Skipping cluster {cluster_id}"

    try:
        if not os.path.exists(cluster_fasta):
            seq_ids = df_cluster['sequence_id'].tolist()
            sequences = load_sequences_with_index(fasta_index, seq_ids) if use_index and fasta_index \
                        else load_sequences_optimized(fasta_path, seq_ids)
            SeqIO.write(sequences, cluster_fasta, "fasta")

        success = align_with_mafft_parallel(cluster_fasta, aligned_fasta, mafft_path, DummyLogger(), threads=2)

        end_time = time.time()
        duration = end_time - start_time
        print(f"Finished cluster {cluster_id}. End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}. Duration: {duration:.2f}s")
        return f"Cluster {cluster_id} - {'aligned' if success else 'failed'} - {duration:.2f}s"
    except Exception as e:
        return f"Error processing cluster {cluster_id}: {str(e)}"


class DummyLogger:
    def write(self, message):
        pass


# --------------------------
# Parallel extraction
# --------------------------
def extract_and_align_clusters_parallel(barcode_path, kmer, mafft_path, logger, max_workers=None, use_index=True):
    output_dir = os.path.join(barcode_path, f"AlignedClusters-k{kmer}")
    os.makedirs(output_dir, exist_ok=True)

    df, clusters = load_top_clusters(barcode_path, kmer)
    if df is None:
        logger.write(f"No cluster CSV found. Exiting.")
        return

    logger.write(f"Found {len(clusters)} clusters: {clusters[:10]}{'...' if len(clusters) > 10 else ''}")

    fasta_file = os.path.join(barcode_path, f"{Path(barcode_path).name}_combined.fasta")
    if not os.path.exists(fasta_file):
        logger.write(f"FASTA file {fasta_file} not found. Exiting.")
        return

    barcode_name = Path(barcode_path).name
    index_path = os.path.join(barcode_path, f"{barcode_name}_fasta_index.pkl")
    fasta_index = create_or_load_fasta_index(fasta_file, index_path, logger) if use_index else None

    process_args = [(cid, df[df['cluster']==cid], fasta_file, output_dir, barcode_name, mafft_path, use_index, fasta_index)
                    for cid in clusters]

    if max_workers is None:
        max_workers = min(mp.cpu_count() // 2, len(clusters), 4)

    logger.write(f"Processing {len(clusters)} clusters using {max_workers} parallel workers...")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(process_cluster, process_args))

    for result in results:
        logger.write(f"{result}")


# --------------------------
# Signal handling for safe exit
# --------------------------
def setup_signal_handlers(logger):
    def handle_exit(signum, frame):
        logger.write(f"Received signal {signum}, flushing logs and exiting...")
        logger.close()
        sys.exit(0)

    signal.signal(signal.SIGTERM, handle_exit)
    signal.signal(signal.SIGINT, handle_exit)


# --------------------------
# Main function
# --------------------------
def main():
    base_path = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus2"
    mafft_path = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Tools\MAFFT\mafft-7.526-win64-signed (1)\mafft-win\mafft.bat"

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    all_barcodes = [f for f in os.listdir(base_path)
                    if os.path.isdir(os.path.join(base_path, f)) and f.lower().startswith('barcode')]
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

    parallel_choice = input("Use parallel processing? (y/n, default=y): ").strip().lower()
    use_parallel = parallel_choice != 'n'
    max_workers = None
    if use_parallel:
        workers_input = input(f"Number of parallel workers (default=auto, max={mp.cpu_count()}): ").strip()
        max_workers = int(workers_input) if workers_input.isdigit() else None

    barcode_path = os.path.join(base_path, barcode_choice)
    log_file_path = os.path.join(barcode_path, f"barcode_clusters_alignment_k{kmer}_{timestamp}.log")
    logger = BufferedLogger(log_file_path)
    setup_signal_handlers(logger)

    try:
        if use_parallel:
            extract_and_align_clusters_parallel(barcode_path, kmer, mafft_path, logger, max_workers)
        else:
            print("Sequential processing not included for brevity")
    finally:
        logger.close()
        print(f"\nLog file saved to: {log_file_path}")


if __name__ == "__main__":
    main()
