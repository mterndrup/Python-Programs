import os
import glob
import datetime
from Bio import SeqIO
import numpy as np
from itertools import product
import umap
import hdbscan
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import psutil
import multiprocessing

# ---------------------------
# Parameters
# ---------------------------
wild2_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Wild4"
log_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Logs"
bases = ['A', 'C', 'G', 'T']
batch_size = 50000  # For batch UMAP processing

os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f"analysis_log_{datetime.datetime.now():%Y%m%d_%H%M%S}.txt")

def log_message(msg):
    """Print and log timestamped message."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"{timestamp} - {msg}"
    print(line, flush=True)
    with open(log_file, "a") as lf:
        lf.write(line + "\n")

def log_memory():
    """Log current memory usage in MB."""
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / (1024 * 1024)
    log_message(f"RAM usage: {mem_mb:.2f} MB")

def kmer_vector(seq, k, kmer_index):
    seq = seq.upper()
    vec = np.zeros(len(kmer_index), dtype=int)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if all(b in bases for b in kmer):
            vec[kmer_index[kmer]] += 1
    return vec

# ---------------------------
# Mode selection
# ---------------------------
log_message("Processing started.")
start_time = datetime.datetime.now()

print("Choose processing mode:")
print("1 = Process each barcode individually")
print("2 = Combine all barcodes into one analysis")
mode = input("Enter 1 or 2: ").strip()

barcode_folders = [d for d in glob.glob(os.path.join(wild2_dir, "*")) if os.path.isdir(d)]

# K-mer selection
if mode == "2":
    while True:
        k = input("Enter k-mer size (6 or 8): ").strip()
        if k in ["6", "8"]:
            k = int(k)
            break
        else:
            print("Invalid choice. Please enter 6 or 8.")
    kmer_list = [k]
else:
    kmer_list = [6, 8]

# ---------------------------
# Process each barcode individually
# ---------------------------
for barcode_dir in barcode_folders:
    barcode_name = os.path.basename(os.path.normpath(barcode_dir))
    raw_dir = os.path.join(barcode_dir, "Raw")
    if not os.path.exists(raw_dir):
        log_message(f"No Raw folder in {barcode_name}, skipping.")
        continue

    fastq_files = glob.glob(os.path.join(raw_dir, "*_trimmed_cutadapt_*"))
    if len(fastq_files) == 0:
        log_message(f"No FASTQ files found in {barcode_name}, skipping.")
        continue

    for k in kmer_list:
        vector_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_vectors.npz")
        if os.path.exists(vector_file):
            log_message(f"{barcode_name} k={k} vectors already exist. Skipping vectorization.")
            continue

        log_message(f"Vectorizing {barcode_name} k={k}...")
        kmers = [''.join(p) for p in product(bases, repeat=k)]
        kmer_index = {kmer: idx for idx, kmer in enumerate(kmers)}

        seq_ids = []
        file_labels = []
        vectors = []

        for fpath in fastq_files:
            log_message(f"Reading file: {os.path.basename(fpath)}")
            for record in SeqIO.parse(fpath, "fastq"):
                vec = kmer_vector(str(record.seq), k, kmer_index)
                vectors.append(vec)
                seq_ids.append(record.id)
                file_labels.append(os.path.basename(fpath))
            log_memory()  # RAM log after each file

        vectors = np.array(vectors, dtype=np.uint16)
        np.savez_compressed(vector_file, vectors=vectors, seq_ids=seq_ids, file_labels=file_labels)
        log_message(f"Saved {vectors.shape[0]} vectors to {vector_file}")
        log_memory()

# ---------------------------
# Combine barcodes (mode 2)
# ---------------------------
if mode == "2":
    log_message("Combined mode selected.")
    combined_plot = os.path.join(log_dir, f"combined_k{k}_UMAP_HDBSCAN.png")
    combined_csv = os.path.join(log_dir, f"combined_k{k}_clusters.csv")
    embeddings_file = os.path.join(log_dir, f"combined_k{k}_embeddings.npy")
    temp_mmap_file = os.path.join(log_dir, f"combined_k{k}_vectors.dat")

    # Determine total vectors and feature size
    total_vectors = 0
    n_features = len([''.join(p) for p in product(bases, repeat=k)])
    for barcode_dir in barcode_folders:
        vector_file = glob.glob(os.path.join(barcode_dir, f"*k{k}_vectors.npz"))
        if vector_file:
            data = np.load(vector_file[0], allow_pickle=True)
            total_vectors += data['vectors'].shape[0]
    log_message(f"Total vectors found: {total_vectors}")
    if total_vectors == 0:
        log_message("No vectors found for combined mode. Exiting.")
        exit()

    # Create or load memory-mapped vectors
    if not os.path.exists(temp_mmap_file):
        log_message(f"Creating memory-mapped file: {temp_mmap_file}")
        combined_vectors = np.memmap(temp_mmap_file, dtype='int32', mode='w+', shape=(total_vectors, n_features))
        all_seq_ids = []
        all_file_labels = []
        idx = 0
        for barcode_dir in barcode_folders:
            vector_file = glob.glob(os.path.join(barcode_dir, f"*k{k}_vectors.npz"))
            if vector_file:
                data = np.load(vector_file[0], allow_pickle=True)
                n = data['vectors'].shape[0]
                combined_vectors[idx:idx + n, :] = data['vectors']
                combined_vectors.flush()
                all_seq_ids.extend(data['seq_ids'])
                all_file_labels.extend(data['file_labels'])
                idx += n
                log_message(f"Loaded {n} vectors from {os.path.basename(vector_file[0])}")
                log_memory()
        combined_vectors.flush()
        log_message("Memory-mapped vectors ready.")
    else:
        log_message(f"Using existing memory-mapped file: {temp_mmap_file}")
        combined_vectors = np.memmap(temp_mmap_file, dtype='int32', mode='r', shape=(total_vectors, n_features))
        all_seq_ids = []
        all_file_labels = []
        for barcode_dir in barcode_folders:
            vector_file = glob.glob(os.path.join(barcode_dir, f"*k{k}_vectors.npz"))
            if vector_file:
                data = np.load(vector_file[0], allow_pickle=True)
                all_seq_ids.extend(data['seq_ids'])
                all_file_labels.extend(data['file_labels'])

    # ---------------------------
    # UMAP: fit once, transform in batches
    # ---------------------------
    if not os.path.exists(embeddings_file):
        log_message("Running UMAP on combined dataset (memory-efficient)...")
        log_memory()
        umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)

        # Fit UMAP on a sample if dataset is huge
        sample_size = min(100000, total_vectors)
        log_message(f"Fitting UMAP on sample of {sample_size} vectors for speed...")
        sample_indices = np.random.choice(total_vectors, sample_size, replace=False)
        umap_model.fit(combined_vectors[sample_indices])

        # Transform full dataset in batches
        embeddings = np.zeros((total_vectors, 2), dtype=np.float32)
        for start in range(0, total_vectors, batch_size):
            end = min(start + batch_size, total_vectors)
            embeddings[start:end] = umap_model.transform(combined_vectors[start:end])
            log_message(f"Transformed vectors {start} to {end}")
            log_memory()

        # Save embeddings
        np.save(embeddings_file, embeddings)
        log_message(f"Saved combined embeddings to {embeddings_file}")
    else:
        log_message(f"Found existing embeddings: {embeddings_file}")
        embeddings = np.load(embeddings_file)

    # ---------------------------
    # HDBSCAN clustering
    # ---------------------------
    if not os.path.exists(combined_csv):
        log_message("Starting HDBSCAN clustering...")
        log_memory()
        clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
        labels = clusterer.fit_predict(embeddings)
        log_memory()
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        log_message(f"Found {n_clusters} clusters.")

        df_clusters = pd.DataFrame({
            'sequence_id': all_seq_ids,
            'source_file': all_file_labels,
            'cluster': labels
        })
        df_clusters.to_csv(combined_csv, index=False)
        log_message(f"Saved combined clusters CSV to {combined_csv}")
    else:
        log_message(f"Found existing CSV: {combined_csv}. Skipping HDBSCAN clustering.")
        df_clusters = pd.read_csv(combined_csv)
        labels = df_clusters['cluster'].values

    # ---------------------------
    # Plot
    # ---------------------------
    log_message("Generating plot...")
    plt.figure(figsize=(12, 10))
    if labels is not None:
        colors = labels
    else:
        colors = np.full(embeddings.shape[0], -1)
    scatter = plt.scatter(embeddings[:, 0], embeddings[:, 1], c=colors, cmap='tab20', s=5, alpha=0.7)
    plt.title(f"UMAP + HDBSCAN Clustering (Combined) (k={k})", fontsize=16)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.colorbar(scatter, label='Cluster ID (-1 = noise)')
    plt.tight_layout()
    plt.savefig(combined_plot, dpi=300)
    plt.close()
    log_message(f"Saved combined plot to {combined_plot}")

# ---------------------------
# Finish
# ---------------------------
end_time = datetime.datetime.now()
total_time = end_time - start_time
log_message(f"Processing finished. Total time: {total_time}")
log_memory()
