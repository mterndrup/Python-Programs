import os
import glob
import datetime
import json
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
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"{timestamp} - {msg}"
    print(line, flush=True)
    with open(log_file, "a") as lf:
        lf.write(line + "\n")


def log_memory():
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / (1024 * 1024)
    log_message(f"RAM usage: {mem_mb:.2f} MB")


def kmer_vector(seq, k, kmer_index):
    seq = seq.upper()
    vec = np.zeros(len(kmer_index), dtype=int)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if all(b in bases for b in kmer):
            vec[kmer_index[kmer]] += 1
    return vec


def save_cluster_colors(barcode_dir, k, colors, seq_ids, labels):
    """Save cluster colors and sequence-to-cluster mapping to JSON file"""
    colors_file = os.path.join(barcode_dir, f"cluster_colors_k{k}.json")

    # Create mapping of sequence IDs to clusters
    seq_cluster_map = {seq_id: int(cluster) for seq_id, cluster in zip(seq_ids, labels)}

    # Convert colors to serializable format (RGBA tuples to lists)
    serializable_colors = {}
    for cluster_id, color in colors.items():
        if isinstance(color, tuple):
            serializable_colors[str(cluster_id)] = list(color)
        else:
            serializable_colors[str(cluster_id)] = color

    color_data = {
        'colors': serializable_colors,
        'sequence_clusters': seq_cluster_map
    }

    with open(colors_file, 'w') as f:
        json.dump(color_data, f, indent=2)

    log_message(f"Saved cluster colors to: {colors_file}")


def load_cluster_colors(barcode_dir, k):
    """Load cluster colors from JSON file"""
    colors_file = os.path.join(barcode_dir, f"cluster_colors_k{k}.json")

    if not os.path.exists(colors_file):
        return None, None

    try:
        with open(colors_file, 'r') as f:
            color_data = json.load(f)

        # Convert string keys back to integers and lists back to tuples
        colors = {}
        for cluster_str, color in color_data['colors'].items():
            cluster_id = int(cluster_str)
            colors[cluster_id] = tuple(color) if isinstance(color, list) else color

        seq_cluster_map = color_data['sequence_clusters']

        log_message(f"Loaded cluster colors from: {colors_file}")
        return colors, seq_cluster_map

    except (json.JSONDecodeError, KeyError) as e:
        log_message(f"Error loading colors from {colors_file}: {e}")
        return None, None


def get_persistent_colors(barcode_dir, k, seq_ids, labels):
    """Get colors for clusters, maintaining consistency across k-mer values"""
    # Try to load colors from k=6 first
    base_colors, base_seq_clusters = load_cluster_colors(barcode_dir, 6)

    if base_colors is None or k == 6:
        # Generate new colors for k=6 or if no base colors exist
        unique_clusters = sorted(c for c in set(labels) if c != -1)
        cmap = plt.get_cmap('tab20')
        colors = {cl: cmap(i % 20) for i, cl in enumerate(unique_clusters)}
        colors[-1] = (0.5, 0.5, 0.5, 0.5)  # grey for noise

        if k == 6:
            # Save colors for k=6 to be used by k=7 and k=8
            save_cluster_colors(barcode_dir, k, colors, seq_ids, labels)

        return colors

    else:
        # For k=7 and k=8, use colors from k=6 where sequences overlap
        current_seq_cluster_map = {seq_id: cluster for seq_id, cluster in zip(seq_ids, labels)}

        # Find which current clusters contain sequences from k=6 clusters
        cluster_mapping = {}
        for seq_id, current_cluster in current_seq_cluster_map.items():
            if seq_id in base_seq_clusters:
                base_cluster = base_seq_clusters[seq_id]
                if current_cluster not in cluster_mapping:
                    cluster_mapping[current_cluster] = {}
                if base_cluster not in cluster_mapping[current_cluster]:
                    cluster_mapping[current_cluster][base_cluster] = 0
                cluster_mapping[current_cluster][base_cluster] += 1

        # Assign colors based on dominant k=6 cluster in each current cluster
        colors = {}
        used_base_colors = set()
        cmap = plt.get_cmap('tab20')
        new_color_idx = len(base_colors)

        unique_clusters = sorted(c for c in set(labels) if c != -1)

        for current_cluster in unique_clusters:
            if current_cluster in cluster_mapping and cluster_mapping[current_cluster]:
                # Find the dominant base cluster
                dominant_base_cluster = max(cluster_mapping[current_cluster],
                                            key=cluster_mapping[current_cluster].get)
                if dominant_base_cluster in base_colors and dominant_base_cluster not in used_base_colors:
                    colors[current_cluster] = base_colors[dominant_base_cluster]
                    used_base_colors.add(dominant_base_cluster)
                else:
                    # Assign new color
                    colors[current_cluster] = cmap(new_color_idx % 20)
                    new_color_idx += 1
            else:
                # New cluster with no k=6 sequences, assign new color
                colors[current_cluster] = cmap(new_color_idx % 20)
                new_color_idx += 1

        colors[-1] = (0.5, 0.5, 0.5, 0.5)  # grey for noise

        return colors


# ---------------------------
# Per-barcode analysis
# ---------------------------
def run_barcode_analysis(k, barcode_dir):
    barcode_name = os.path.basename(os.path.normpath(barcode_dir))
    vector_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_vectors.npz")
    clusters_csv = os.path.join(barcode_dir, f"{barcode_name}_k{k}_clusters.csv")
    embeddings_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_UMAP.npy")
    plot_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_UMAP_HDBSCAN.png")

    if not os.path.exists(vector_file):
        log_message(f"Vectors not found for {barcode_name} k={k}, skipping analysis.")
        return

    log_message(f"Loading vectors for {barcode_name} k={k}...")
    data = np.load(vector_file, allow_pickle=True)
    vectors = data['vectors']
    seq_ids = data['seq_ids']
    file_labels = data['file_labels']

    # ---------------------------
    # UMAP
    # ---------------------------
    if not os.path.exists(embeddings_file):
        log_message(f"Running UMAP for {barcode_name} k={k}...")
        umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
        embeddings = umap_model.fit_transform(vectors)
        np.save(embeddings_file, embeddings)
        log_message(f"Saved UMAP embeddings: {embeddings_file}")
    else:
        log_message(f"Found existing UMAP embeddings: {embeddings_file}")
        embeddings = np.load(embeddings_file)

    # ---------------------------
    # HDBSCAN
    # ---------------------------
    if not os.path.exists(clusters_csv):
        log_message(f"Clustering {barcode_name} with HDBSCAN...")
        clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
        labels = clusterer.fit_predict(embeddings)
        df_clusters = pd.DataFrame({
            'sequence_id': seq_ids,
            'source_file': file_labels,
            'cluster': labels
        })
        df_clusters.to_csv(clusters_csv, index=False)
        log_message(f"Saved clusters CSV: {clusters_csv}")
    else:
        log_message(f"Found existing clusters CSV: {clusters_csv}")
        df_clusters = pd.read_csv(clusters_csv)
        labels = df_clusters['cluster'].values

    # ---------------------------
    # Get persistent colors
    # ---------------------------
    colors = get_persistent_colors(barcode_dir, k, seq_ids, labels)
    point_colors = [colors.get(cl, (0.8, 0.8, 0.8, 0.5)) for cl in labels]

    # ---------------------------
    # Plot
    # ---------------------------
    plt.figure(figsize=(12, 10))
    plt.scatter(embeddings[:, 0], embeddings[:, 1], c=point_colors, s=5, alpha=0.7)
    plt.title(f"UMAP + HDBSCAN ({barcode_name}, k={k})", fontsize=16)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()
    log_message(f"Saved plot: {plot_file}")


# ---------------------------
# Combined analysis
# ---------------------------
def run_combined_analysis(k, barcode_dirs, base_colors=None):
    log_message(f"Combined mode analysis for k={k}")
    combined_plot = os.path.join(log_dir, f"combined_k{k}_UMAP_HDBSCAN.png")
    combined_csv = os.path.join(log_dir, f"combined_k{k}_clusters.csv")
    embeddings_file = os.path.join(log_dir, f"combined_k{k}_embeddings.npy")
    temp_mmap_file = os.path.join(log_dir, f"combined_k{k}_vectors.dat")

    n_features = len([''.join(p) for p in product(bases, repeat=k)])
    total_vectors = 0
    for barcode_dir in barcode_dirs:
        vector_file = glob.glob(os.path.join(barcode_dir, f"*k{k}_vectors.npz"))
        if vector_file:
            data = np.load(vector_file[0], allow_pickle=True)
            total_vectors += data['vectors'].shape[0]

    log_message(f"Total vectors found for k={k}: {total_vectors}")
    if total_vectors == 0:
        log_message("No vectors found. Exiting.")
        return None

    # ---------------------------
    # Load/Create memmap
    # ---------------------------
    if not os.path.exists(temp_mmap_file):
        log_message(f"Creating memmap file: {temp_mmap_file}")
        combined_vectors = np.memmap(temp_mmap_file, dtype='int32', mode='w+', shape=(total_vectors, n_features))
        all_seq_ids, all_file_labels = [], []
        idx = 0
        for barcode_dir in barcode_dirs:
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
    else:
        log_message(f"Using existing memmap: {temp_mmap_file}")
        combined_vectors = np.memmap(temp_mmap_file, dtype='int32', mode='r', shape=(total_vectors, n_features))
        all_seq_ids, all_file_labels = [], []
        for barcode_dir in barcode_dirs:
            vector_file = glob.glob(os.path.join(barcode_dir, f"*k{k}_vectors.npz"))
            if vector_file:
                data = np.load(vector_file[0], allow_pickle=True)
                all_seq_ids.extend(data['seq_ids'])
                all_file_labels.extend(data['file_labels'])

    # ---------------------------
    # UMAP
    # ---------------------------
    if not os.path.exists(embeddings_file):
        log_message("Running UMAP (memory-efficient)...")
        log_memory()
        umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
        sample_size = min(100000, total_vectors)
        log_message(f"Fitting UMAP on {sample_size} sample points...")
        sample_indices = np.random.choice(total_vectors, sample_size, replace=False)
        umap_model.fit(combined_vectors[sample_indices])
        embeddings = np.zeros((total_vectors, 2), dtype=np.float32)
        for start in range(0, total_vectors, batch_size):
            end = min(start + batch_size, total_vectors)
            embeddings[start:end] = umap_model.transform(combined_vectors[start:end])
            log_message(f"Transformed vectors {start} to {end}")
            log_memory()
        np.save(embeddings_file, embeddings)
        log_message(f"Saved embeddings: {embeddings_file}")
    else:
        log_message(f"Found existing embeddings: {embeddings_file}")
        embeddings = np.load(embeddings_file)

    # ---------------------------
    # HDBSCAN
    # ---------------------------
    if not os.path.exists(combined_csv):
        log_message("Clustering with HDBSCAN...")
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
        log_message(f"Saved clusters CSV: {combined_csv}")
    else:
        log_message(f"Found existing CSV: {combined_csv}")
        df_clusters = pd.read_csv(combined_csv)
        labels = df_clusters['cluster'].values

    # ---------------------------
    # Plotting with persistent colors
    # ---------------------------
    if base_colors is None:
        base_colors = {}

    unique_clusters = sorted(c for c in set(labels) if c != -1)
    cmap = plt.get_cmap('tab20')
    next_color_idx = len(base_colors)

    for cl in unique_clusters:
        if cl not in base_colors:
            base_colors[cl] = cmap(next_color_idx % 20)
            next_color_idx += 1
    base_colors[-1] = (0.5, 0.5, 0.5, 0.5)  # grey for noise

    colors = {cl: base_colors.get(cl, (0.8, 0.8, 0.8, 0.5)) for cl in set(labels)}

    log_message(f"Generating plot for k={k}...")
    plt.figure(figsize=(12, 10))
    point_colors = [colors.get(cl, (0.8, 0.8, 0.8, 0.5)) for cl in labels]
    plt.scatter(embeddings[:, 0], embeddings[:, 1], c=point_colors, s=5, alpha=0.7)
    plt.title(f"UMAP + HDBSCAN (Combined) (k={k})", fontsize=16)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(combined_plot, dpi=300)
    plt.close()
    log_message(f"Saved plot: {combined_plot}")

    return base_colors


# ---------------------------
# Start
# ---------------------------
log_message("Processing started.")
start_time = datetime.datetime.now()

print("Choose processing mode:")
print("1 = Process each barcode individually")
print("2 = Combine all barcodes into one analysis")
mode = input("Enter 1 or 2: ").strip()

barcode_folders = [d for d in glob.glob(os.path.join(wild2_dir, "*")) if os.path.isdir(d)]

# ---------------------------
# Mode 1: vectorization & per-barcode analysis
# ---------------------------
if mode == "1":
    kmer_list = [6, 7, 8]
    base_colors = None
    for k in kmer_list:
        for barcode_dir in barcode_folders:
            barcode_name = os.path.basename(os.path.normpath(barcode_dir))
            raw_dir = os.path.join(barcode_dir, "Raw")
            if not os.path.exists(raw_dir):
                log_message(f"No Raw folder in {barcode_name}, skipping.")
                continue

            fastq_files = glob.glob(os.path.join(raw_dir, "*_trimmed_cutadapt_*"))
            if not fastq_files:
                log_message(f"No FASTQ files in {barcode_name}, skipping.")
                continue

            vector_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_vectors.npz")
            if os.path.exists(vector_file):
                log_message(f"{barcode_name} k={k} vectors exist. Skipping vectorization.")
            else:
                log_message(f"Vectorizing {barcode_name} k={k}...")
                kmers = [''.join(p) for p in product(bases, repeat=k)]
                kmer_index = {km: idx for idx, km in enumerate(kmers)}

                seq_ids, file_labels, vectors = [], [], []
                for fpath in fastq_files:
                    log_message(f"Reading {os.path.basename(fpath)}")
                    for record in SeqIO.parse(fpath, "fastq"):
                        vec = kmer_vector(str(record.seq), k, kmer_index)
                        vectors.append(vec)
                        seq_ids.append(record.id)
                        file_labels.append(os.path.basename(fpath))
                    log_memory()

                vectors = np.array(vectors, dtype=np.uint16)
                np.savez_compressed(vector_file, vectors=vectors, seq_ids=seq_ids, file_labels=file_labels)
                log_message(f"Saved {vectors.shape[0]} vectors to {vector_file}")

            # Run per-barcode analysis
            run_barcode_analysis(k, barcode_dir)

# ---------------------------
# Mode 2: combined analysis
# ---------------------------
elif mode == "2":
    kmer_list = [6, 7, 8]
    base_colors = None
    for k in kmer_list:
        base_colors = run_combined_analysis(k, barcode_folders, base_colors)

end_time = datetime.datetime.now()
log_message(f"Processing finished. Total time: {end_time - start_time}")