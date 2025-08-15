import os
import glob
import datetime
from Bio import SeqIO
import numpy as np
from itertools import product
import pandas as pd
import psutil

# ---------------------------
# Parameters
# ---------------------------
wild2_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Wild2"
log_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Logs"
bases = ['A', 'C', 'G', 'T']
batch_size = 50000  # For batch processing in combined mode

os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f"analysis_log_{datetime.datetime.now():%Y%m%d_%H%M%S}.txt")

# ---------------------------
# Logging Functions
# ---------------------------
def log_message(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"{timestamp} - {msg}"
    print(line, flush=True)
    with open(log_file, "a") as lf:
        lf.write(line + "\n")

def log_memory():
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / (1024 * 1024)
    log_message(f"RAM usage: {mem_mb:.2f} MB")

# ---------------------------
# K-mer Vectorization
# ---------------------------
def kmer_vector(seq, k, kmer_index):
    seq = seq.upper()
    vec = np.zeros(len(kmer_index), dtype=int)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if all(b in bases for b in kmer):
            vec[kmer_index[kmer]] += 1
    return vec

# ---------------------------
# User Mode Selection
# ---------------------------
log_message("Processing started.")
print("Choose processing mode:")
print("1 = Process all barcodes individually")
print("2 = Combine all barcodes into one analysis")
print("3 = Process individually selected barcodes")
mode = input("Enter 1, 2, or 3: ").strip()

barcode_folders = sorted([d for d in glob.glob(os.path.join(wild2_dir, "*")) if os.path.isdir(d)])

if mode == "3":
    print("Available barcodes:")
    for folder in barcode_folders:
        print(os.path.basename(folder))
    selected = input("Enter comma-separated barcode names to process (e.g., 77,78): ").strip()
    selected_names = [x.strip() for x in selected.split(",") if x.strip()]
    barcode_folders = [folder for folder in barcode_folders if os.path.basename(folder) in selected_names]

# K-mer selection
if mode == "2":
    while True:
        k = input("Enter k-mer size (6, 7, or 8): ").strip()
        if k in ["6", "7", "8"]:
            kmer_list = [int(k)]
            break
        else:
            print("Invalid choice. Please enter 6, 7, or 8.")
else:
    kmer_list = [6, 7, 8]

# ---------------------------
# Import Heavy Libraries
# ---------------------------
import umap
import hdbscan
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------
# Function: Process Single Barcode
# ---------------------------
def process_barcode(barcode_dir, kmer_list):
    barcode_name = os.path.basename(barcode_dir)
    raw_dir = os.path.join(barcode_dir, "Raw")
    if not os.path.exists(raw_dir):
        log_message(f"No Raw folder in {barcode_name}, skipping.")
        return

    fastq_files = glob.glob(os.path.join(raw_dir, "*_trimmed_cutadapt_*"))
    if len(fastq_files) == 0:
        log_message(f"No FASTQ files found in {barcode_name}.")
        return

    charts_dir = os.path.join(barcode_dir, "Charts")
    os.makedirs(charts_dir, exist_ok=True)

    # k6 color file path
    k6_color_file = os.path.join(barcode_dir, f"k6_label_colors.npy")
    k6_colors_dict = None
    if os.path.exists(k6_color_file):
        k6_colors_dict = np.load(k6_color_file, allow_pickle=True).item()

    for k in kmer_list:
        vector_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_vectors.npz")
        umap_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_UMAP.npy")
        cluster_csv = os.path.join(barcode_dir, f"{barcode_name}_k{k}_clusters.csv")
        plot_file = os.path.join(charts_dir, f"{barcode_name}_k{k}_UMAP_HDBSCAN.png")

        # ---------------------------
        # Vectorization
        # ---------------------------
        if os.path.exists(vector_file):
            log_message(f"Vectors file exists for {barcode_name} k={k}, skipping FASTQ reading: {vector_file}")
            data = np.load(vector_file, allow_pickle=True)
            vectors = data['vectors']
            seq_ids = data['seq_ids']
            file_labels = data['file_labels']
        else:
            log_message(f"Vector file missing, processing FASTQ for {barcode_name} k={k}...")
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
                log_memory()

            vectors = np.array(vectors, dtype=np.uint16)
            np.savez_compressed(vector_file, vectors=vectors, seq_ids=seq_ids, file_labels=file_labels)
            log_message(f"Saved {vectors.shape[0]} vectors to {vector_file}")
            log_memory()

        # ---------------------------
        # UMAP embeddings
        # ---------------------------
        if os.path.exists(umap_file):
            log_message(f"Loading existing UMAP embeddings: {umap_file}")
            embeddings = np.load(umap_file)
        else:
            log_message(f"Running UMAP for {barcode_name} k={k}...")
            log_memory()
            umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
            embeddings = umap_model.fit_transform(vectors)
            np.save(umap_file, embeddings)
            log_message(f"Saved UMAP embeddings to {umap_file}")
            log_memory()

        # ---------------------------
        # HDBSCAN clustering
        # ---------------------------
        if os.path.exists(cluster_csv):
            log_message(f"Loading existing cluster CSV: {cluster_csv}")
            df_clusters = pd.read_csv(cluster_csv)
            labels = df_clusters['cluster'].values
        else:
            log_message(f"Running HDBSCAN clustering for {barcode_name} k={k}...")
            log_memory()
            clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
            labels = clusterer.fit_predict(embeddings)
            df_clusters = pd.DataFrame({
                'sequence_id': seq_ids,
                'source_file': file_labels,
                'cluster': labels
            })
            df_clusters.to_csv(cluster_csv, index=False)
            log_message(f"Saved clusters CSV to {cluster_csv}")
            log_memory()

        # ---------------------------
        # Plot generation with retroactive k6 colors
        # ---------------------------
        if k == 6 and k6_colors_dict is None:
            unique_clusters = np.unique(labels)
            rng = np.random.default_rng(42)
            colors = rng.random((len(unique_clusters), 4))
            colors[:, 3] = 1  # full alpha
            k6_colors_dict = dict(zip(unique_clusters, colors))
            np.save(k6_color_file, k6_colors_dict)
            log_message(f"Saved k6 color mapping to {k6_color_file}")

        if k != 6 and k6_colors_dict is not None:
            # Use k6 colors where cluster label matches; else gray
            cluster_colors = np.array([k6_colors_dict.get(lbl, np.array([0.5,0.5,0.5,1])) for lbl in labels])
        else:
            # fallback random colors
            unique_clusters = np.unique(labels)
            rng = np.random.default_rng(42)
            colors = rng.random((len(unique_clusters), 4))
            colors[:, 3] = 1
            cluster_colors = np.array([dict(zip(unique_clusters, colors)).get(lbl) for lbl in labels])

        plt.figure(figsize=(12, 10))
        plt.scatter(embeddings[:, 0], embeddings[:, 1], c=cluster_colors, s=5, alpha=0.7)
        plt.title(f"{barcode_name} UMAP + HDBSCAN (k={k})", fontsize=16)
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.colorbar(label='Cluster ID (-1 = noise)')
        plt.tight_layout()
        plt.savefig(plot_file, dpi=300)
        plt.close()
        log_message(f"Saved plot to {plot_file}")
        log_memory()

# ---------------------------
# Process Individual or Selected Barcodes
# ---------------------------
if mode in ["1", "3"]:
    for barcode_dir in barcode_folders:
        process_barcode(barcode_dir, kmer_list)

# ---------------------------
# Combined Mode
# ---------------------------
if mode == "2":
    log_message("Combined mode selected.")
    timestamp_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    charts_dir = os.path.join(r"C:\Users\ketgl\OneDrive\Desktop\Sandbox", "Charts")
    os.makedirs(charts_dir, exist_ok=True)
    combined_plot = os.path.join(charts_dir, f"combined_k{k}_UMAP_HDBSCAN_{timestamp_str}.png")
    combined_csv = os.path.join(log_dir, f"combined_k{k}_clusters_{timestamp_str}.csv")
    embeddings_file = os.path.join(log_dir, f"combined_k{k}_embeddings_{timestamp_str}.npy")
    temp_mmap_file = os.path.join(log_dir, f"combined_k{k}_vectors_{timestamp_str}.dat")

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

    log_message(f"Creating memory-mapped file: {temp_mmap_file}")
    combined_vectors = np.memmap(temp_mmap_file, dtype='int32', mode='w+', shape=(total_vectors, n_features))
    all_seq_ids = []
    all_file_labels = []
    all_labels_k6 = []
    idx = 0

    # Collect k6 colors for combined mode if available
    k6_colors_combined = None
    for barcode_dir in barcode_folders:
        k6_file = os.path.join(barcode_dir, "k6_label_colors.npy")
        if os.path.exists(k6_file):
            k6_colors_combined = np.load(k6_file, allow_pickle=True).item()
            break

    for barcode_dir in barcode_folders:
        vector_file = glob.glob(os.path.join(barcode_dir, f"*k{k}_vectors.npz"))
        if not vector_file:
            continue
        data = np.load(vector_file[0], allow_pickle=True)
        vectors = data['vectors']
        n = vectors.shape[0]
        combined_vectors[idx:idx+n, :] = vectors
        all_seq_ids.extend(data['seq_ids'])
        all_file_labels.extend(data['file_labels'])
        idx += n
    log_memory()

    # UMAP
    log_message("Running combined UMAP...")
    umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
    combined_embeddings = umap_model.fit_transform(combined_vectors)
    np.save(embeddings_file, combined_embeddings)
    log_message(f"Saved combined embeddings to {embeddings_file}")
    log_memory()

    # HDBSCAN
    log_message("Running combined HDBSCAN...")
    clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
    labels = clusterer.fit_predict(combined_embeddings)
    df_combined = pd.DataFrame({
        'sequence_id': all_seq_ids,
        'source_file': all_file_labels,
        'cluster': labels
    })
    df_combined.to_csv(combined_csv, index=False)
    log_message(f"Saved combined clusters CSV to {combined_csv}")
    log_memory()

    # Plot
    plt.figure(figsize=(12, 10))
    if k6_colors_combined is not None:
        cluster_colors = np.array([k6_colors_combined.get(lbl, np.array([0.5,0.5,0.5,1])) for lbl in labels])
    else:
        unique_clusters = np.unique(labels)
        rng = np.random.default_rng(42)
        colors = rng.random((len(unique_clusters), 4))
        colors[:, 3] = 1
        color_dict = dict(zip(unique_clusters, colors))
        cluster_colors = np.array([color_dict.get(lbl, np.array([0.5,0.5,0.5,1])) for lbl in labels])

    plt.scatter(combined_embeddings[:, 0], combined_embeddings[:, 1], c=cluster_colors, s=5, alpha=0.7)
    plt.title(f"Combined UMAP + HDBSCAN (k={k})", fontsize=16)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.colorbar(label='Cluster ID (-1 = noise)')
    plt.tight_layout()
    plt.savefig(combined_plot, dpi=300)
    plt.close()
    log_message(f"Saved combined plot to {combined_plot}")
    log_memory()

log_message("Processing complete.")
