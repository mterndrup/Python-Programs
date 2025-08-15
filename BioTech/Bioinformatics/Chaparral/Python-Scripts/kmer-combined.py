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
import gc
from collections import defaultdict

# ---------------------------
# Parameters
# ---------------------------
wild2_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Wild4"
log_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Logs"
charts_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Charts"
bases = ['A', 'C', 'G', 'T']
batch_size = 10000  # Reduced batch size for memory efficiency
memory_limit_mb = 8000  # Set memory limit (adjust based on your system)

os.makedirs(log_dir, exist_ok=True)
os.makedirs(charts_dir, exist_ok=True)
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
    return mem_mb


def check_memory_limit():
    """Check if memory usage is approaching limit"""
    mem_mb = log_memory()
    if mem_mb > memory_limit_mb * 0.9:  # 90% of limit
        log_message(f"WARNING: Memory usage ({mem_mb:.2f} MB) approaching limit ({memory_limit_mb} MB)")
        return True
    return False


def force_garbage_collection():
    """Force garbage collection and log memory reduction"""
    mem_before = log_memory()
    gc.collect()
    mem_after = log_memory()
    log_message(f"Garbage collection: {mem_before:.2f} -> {mem_after:.2f} MB (freed {mem_before - mem_after:.2f} MB)")


def sparse_kmer_vector(seq, k, min_frequency=2):
    """Create sparse k-mer vector using dictionary, filtering low-frequency k-mers"""
    seq = seq.upper()
    kmer_counts = defaultdict(int)

    # Count k-mers
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if all(b in bases for b in kmer):
            kmer_counts[kmer] += 1

    # Filter low-frequency k-mers to reduce dimensionality
    filtered_counts = {kmer: count for kmer, count in kmer_counts.items()
                       if count >= min_frequency}

    return filtered_counts


def create_sparse_vector_batch(sequences, k, kmer_index=None, min_frequency=2):
    """Process sequences in batches to create sparse vectors"""
    if kmer_index is None:
        # First pass: collect all k-mers across sequences to build index
        all_kmers = set()
        for seq in sequences:
            sparse_counts = sparse_kmer_vector(seq, k, min_frequency)
            all_kmers.update(sparse_counts.keys())

        kmer_index = {kmer: idx for idx, kmer in enumerate(sorted(all_kmers))}
        log_message(f"Created k-mer index with {len(kmer_index)} unique k-mers (filtered, min_freq={min_frequency})")

    # Second pass: create vectors
    vectors = []
    for seq in sequences:
        sparse_counts = sparse_kmer_vector(seq, k, min_frequency)
        vec = np.zeros(len(kmer_index), dtype=np.uint16)
        for kmer, count in sparse_counts.items():
            if kmer in kmer_index:
                vec[kmer_index[kmer]] = min(count, 65535)  # Cap at uint16 max
        vectors.append(vec)

    return np.array(vectors), kmer_index


def process_sequences_streaming(fastq_files, k, vector_file, chunk_size=1000):
    """Process sequences in streaming fashion to reduce memory usage"""
    log_message(f"Streaming processing with chunk_size={chunk_size}")

    # First pass: collect k-mer statistics
    log_message("First pass: collecting k-mer statistics...")
    kmer_freq = defaultdict(int)
    total_seqs = 0

    for fpath in fastq_files:
        log_message(f"Scanning {os.path.basename(fpath)} for k-mer frequencies...")
        sequences_chunk = []

        for record in SeqIO.parse(fpath, "fastq"):
            sequences_chunk.append(str(record.seq))
            total_seqs += 1

            if len(sequences_chunk) >= chunk_size:
                # Process chunk
                for seq in sequences_chunk:
                    sparse_counts = sparse_kmer_vector(seq, k, min_frequency=1)
                    for kmer, count in sparse_counts.items():
                        kmer_freq[kmer] += count

                sequences_chunk = []

                if total_seqs % (chunk_size * 10) == 0:
                    log_memory()

        # Process remaining sequences
        if sequences_chunk:
            for seq in sequences_chunk:
                sparse_counts = sparse_kmer_vector(seq, k, min_frequency=1)
                for kmer, count in sparse_counts.items():
                    kmer_freq[kmer] += count

    log_message(f"Found {len(kmer_freq)} unique k-mers across {total_seqs} sequences")

    # Filter k-mers by frequency (remove very rare ones to reduce dimensionality)
    min_global_freq = max(2, total_seqs // 10000)  # Adaptive threshold
    filtered_kmers = {kmer for kmer, freq in kmer_freq.items() if freq >= min_global_freq}
    kmer_index = {kmer: idx for idx, kmer in enumerate(sorted(filtered_kmers))}

    log_message(f"Filtered to {len(kmer_index)} k-mers (min_global_freq={min_global_freq})")

    # Clear memory
    del kmer_freq
    force_garbage_collection()

    # Second pass: create vectors
    log_message("Second pass: creating vectors...")
    all_vectors = []
    all_seq_ids = []
    all_file_labels = []
    processed_seqs = 0

    for fpath in fastq_files:
        log_message(f"Vectorizing {os.path.basename(fpath)}...")
        sequences_chunk = []
        seq_ids_chunk = []
        file_labels_chunk = []

        for record in SeqIO.parse(fpath, "fastq"):
            sequences_chunk.append(str(record.seq))
            seq_ids_chunk.append(record.id)
            file_labels_chunk.append(os.path.basename(fpath))
            processed_seqs += 1

            if len(sequences_chunk) >= chunk_size:
                # Process chunk
                vectors_chunk, _ = create_sparse_vector_batch(sequences_chunk, k, kmer_index)
                all_vectors.append(vectors_chunk)
                all_seq_ids.extend(seq_ids_chunk)
                all_file_labels.extend(file_labels_chunk)

                # Clear chunk data
                sequences_chunk = []
                seq_ids_chunk = []
                file_labels_chunk = []

                if processed_seqs % (chunk_size * 5) == 0:
                    if check_memory_limit():
                        force_garbage_collection()

                log_message(f"Processed {processed_seqs}/{total_seqs} sequences")

        # Process remaining sequences
        if sequences_chunk:
            vectors_chunk, _ = create_sparse_vector_batch(sequences_chunk, k, kmer_index)
            all_vectors.append(vectors_chunk)
            all_seq_ids.extend(seq_ids_chunk)
            all_file_labels.extend(file_labels_chunk)

    # Combine all vectors
    log_message("Combining all vectors...")
    if all_vectors:
        combined_vectors = np.vstack(all_vectors)
        log_message(f"Final vector shape: {combined_vectors.shape}")

        # Save to compressed format
        np.savez_compressed(vector_file,
                            vectors=combined_vectors,
                            seq_ids=all_seq_ids,
                            file_labels=all_file_labels,
                            kmer_index=kmer_index)

        log_message(f"Saved {combined_vectors.shape[0]} vectors to {vector_file}")
        return combined_vectors.shape[0]
    else:
        log_message("No vectors created")
        return 0


def incremental_umap(vectors, embeddings_file, sample_ratio=0.1, batch_size=5000):
    """Run UMAP with incremental processing for large datasets"""
    n_samples = vectors.shape[0]

    if n_samples > 50000:  # Use incremental processing for large datasets
        log_message(f"Large dataset ({n_samples} samples), using incremental UMAP...")

        # Train on subsample
        sample_size = max(int(n_samples * sample_ratio), 10000)
        sample_size = min(sample_size, 50000)  # Cap at 50k

        log_message(f"Training UMAP on {sample_size} samples...")
        sample_indices = np.random.choice(n_samples, sample_size, replace=False)

        umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2,
                               random_state=42, low_memory=True)

        if hasattr(vectors, 'shape'):
            sample_data = vectors[sample_indices]
        else:
            sample_data = vectors  # Already a sample

        umap_model.fit(sample_data)

        # Transform in batches
        embeddings = np.zeros((n_samples, 2), dtype=np.float32)

        for start in range(0, n_samples, batch_size):
            end = min(start + batch_size, n_samples)

            if hasattr(vectors, 'shape'):
                batch_data = vectors[start:end]
            else:
                batch_data = vectors[start:end]

            embeddings[start:end] = umap_model.transform(batch_data)

            log_message(f"UMAP transform: {start}-{end}/{n_samples}")

            if check_memory_limit():
                force_garbage_collection()

        np.save(embeddings_file, embeddings)
        log_message(f"Saved incremental UMAP embeddings: {embeddings_file}")

    else:
        # Standard UMAP for smaller datasets
        log_message(f"Standard UMAP for {n_samples} samples...")
        umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2,
                               random_state=42, low_memory=True)
        embeddings = umap_model.fit_transform(vectors)
        np.save(embeddings_file, embeddings)
        log_message(f"Saved UMAP embeddings: {embeddings_file}")

    return embeddings


def memory_efficient_hdbscan(embeddings, min_cluster_size=10):
    """Run HDBSCAN with memory optimization"""
    log_message("Running memory-efficient HDBSCAN...")

    # For very large datasets, consider using approximate clustering
    if embeddings.shape[0] > 100000:
        log_message("Large dataset detected, using approximate clustering...")
        # Sample for initial clustering
        sample_size = 50000
        sample_indices = np.random.choice(embeddings.shape[0], sample_size, replace=False)
        sample_embeddings = embeddings[sample_indices]

        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size,
                                    approx_min_span_tree=True,
                                    gen_min_span_tree=False)
        sample_labels = clusterer.fit_predict(sample_embeddings)

        # Predict labels for remaining points
        # This is a simplified approach - in practice, you might want more sophisticated methods
        from sklearn.neighbors import NearestNeighbors

        labeled_points = sample_embeddings[sample_labels != -1]
        labeled_labels = sample_labels[sample_labels != -1]

        if len(labeled_points) > 0:
            nn = NearestNeighbors(n_neighbors=1)
            nn.fit(labeled_points)

            all_labels = np.full(embeddings.shape[0], -1)
            all_labels[sample_indices] = sample_labels

            unlabeled_mask = all_labels == -1
            if np.any(unlabeled_mask):
                distances, indices = nn.kneighbors(embeddings[unlabeled_mask])
                # Assign labels based on nearest labeled neighbor (with distance threshold)
                threshold = np.percentile(distances, 50)  # Use median distance as threshold
                close_enough = distances.flatten() < threshold
                all_labels[np.where(unlabeled_mask)[0][close_enough]] = labeled_labels[indices.flatten()[close_enough]]

            labels = all_labels
        else:
            labels = np.full(embeddings.shape[0], -1)

    else:
        # Standard HDBSCAN for smaller datasets
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
        labels = clusterer.fit_predict(embeddings)

    return labels


# ---------------------------
# Updated per-barcode analysis
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

    log_message(f"Loaded {vectors.shape[0]} vectors with {vectors.shape[1]} features")

    # ---------------------------
    # UMAP with memory efficiency
    # ---------------------------
    if not os.path.exists(embeddings_file):
        embeddings = incremental_umap(vectors, embeddings_file)
    else:
        log_message(f"Found existing UMAP embeddings: {embeddings_file}")
        embeddings = np.load(embeddings_file)

    # Clear vectors from memory if not needed
    del vectors
    force_garbage_collection()

    # ---------------------------
    # HDBSCAN with memory efficiency
    # ---------------------------
    if not os.path.exists(clusters_csv):
        labels = memory_efficient_hdbscan(embeddings)

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
    # Get persistent colors (unchanged)
    # ---------------------------
    colors = get_persistent_colors(barcode_dir, k, seq_ids, labels)
    point_colors = [colors.get(cl, (0.8, 0.8, 0.8, 0.5)) for cl in labels]

    # ---------------------------
    # Plot (unchanged)
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


# [Include the rest of your original functions with similar memory optimizations...]
# [The save_cluster_colors, load_cluster_colors, get_persistent_colors functions remain the same]

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
# Start (updated with memory-efficient vectorization)
# ---------------------------
log_message("Processing started.")
start_time = datetime.datetime.now()

print("Choose processing mode:")
print("1 = Process each barcode individually")
print("2 = Combine all barcodes into one analysis")
mode = input("Enter 1 or 2: ").strip()

barcode_folders = [d for d in glob.glob(os.path.join(wild2_dir, "*")) if os.path.isdir(d)]

# ---------------------------
# Mode 1: memory-efficient vectorization & per-barcode analysis
# ---------------------------
if mode == "1":
    # Ask if this is a restart
    print("\nIs this a restart from a previous run?")
    print("1 = Start from beginning")
    print("2 = Restart from specific barcode and k-mer")
    restart_mode = input("Enter 1 or 2: ").strip()

    start_barcode_idx = 0
    start_k_idx = 0
    kmer_list = [6, 7, 8]

    if restart_mode == "2":
        # Show available barcodes
        print("\nAvailable barcodes:")
        barcode_names = [os.path.basename(os.path.normpath(d)) for d in barcode_folders]
        for i, name in enumerate(barcode_names):
            print(f"{i + 1}. {name}")

        # Get starting barcode
        while True:
            try:
                barcode_input = input(f"\nEnter barcode number to start from (1-{len(barcode_names)}): ").strip()
                start_barcode_idx = int(barcode_input) - 1
                if 0 <= start_barcode_idx < len(barcode_names):
                    break
                else:
                    print(f"Please enter a number between 1 and {len(barcode_names)}")
            except ValueError:
                print("Please enter a valid number")

        # Get starting k-mer value
        print("\nAvailable k-mer values:")
        for i, k in enumerate(kmer_list):
            print(f"{i + 1}. k={k}")

        while True:
            try:
                k_input = input(f"Enter k-mer number to start from (1-{len(kmer_list)}): ").strip()
                start_k_idx = int(k_input) - 1
                if 0 <= start_k_idx < len(kmer_list):
                    break
                else:
                    print(f"Please enter a number between 1 and {len(kmer_list)}")
            except ValueError:
                print("Please enter a valid number")

        selected_barcode = barcode_names[start_barcode_idx]
        selected_k = kmer_list[start_k_idx]

        print(f"\nRestarting from barcode '{selected_barcode}' with k={selected_k}")
        log_message(
            f"RESTART: Starting from barcode {selected_barcode} (index {start_barcode_idx}) with k={selected_k} (index {start_k_idx})")

        # Confirm before proceeding
        confirm = input("Continue with this selection? (y/n): ").strip().lower()
        if confirm != 'y':
            print("Restart cancelled.")
            exit()

    # Main processing loop with restart capability
    for k_idx, k in enumerate(kmer_list):
        if k_idx < start_k_idx:
            log_message(f"Skipping k={k} (before restart point)")
            continue

        for barcode_idx, barcode_dir in enumerate(barcode_folders):
            # Skip barcodes before restart point for the current k-mer
            if k_idx == start_k_idx and barcode_idx < start_barcode_idx:
                continue
            elif k_idx > start_k_idx:
                # For k-mer values after restart point, process all barcodes
                pass

            barcode_name = os.path.basename(os.path.normpath(barcode_dir))
            log_message(f"Processing barcode {barcode_idx + 1}/{len(barcode_folders)}: {barcode_name} with k={k}")

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
                log_message(f"Memory-efficient vectorization for {barcode_name} k={k}...")

                try:
                    # Use streaming processing for memory efficiency
                    num_vectors = process_sequences_streaming(fastq_files, k, vector_file)

                    if num_vectors == 0:
                        log_message(f"No vectors created for {barcode_name} k={k}")
                        continue

                    force_garbage_collection()

                except Exception as e:
                    log_message(f"Error during vectorization for {barcode_name} k={k}: {str(e)}")
                    continue

            # Run per-barcode analysis with memory optimizations
            try:
                run_barcode_analysis(k, barcode_dir)
                force_garbage_collection()
                log_message(f"Successfully completed {barcode_name} k={k}")

            except MemoryError as e:
                log_message(f"Memory error processing {barcode_name} k={k}: {str(e)}")
                log_message("Consider reducing chunk_size or memory_limit_mb")
                continue
            except Exception as e:
                log_message(f"Unexpected error processing {barcode_name} k={k}: {str(e)}")
                continue

# ... (rest of the script for mode 2 would follow similar memory optimization patterns)

end_time = datetime.datetime.now()
log_message(f"Processing finished. Total time: {end_time - start_time}")