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


def adaptive_min_cluster_size(n_samples, base_min_size=10):
    """Calculate adaptive minimum cluster size based on dataset size"""
    if n_samples < 1000:
        return base_min_size
    elif n_samples < 10000:
        return max(base_min_size, n_samples // 200)
    elif n_samples < 100000:
        return max(base_min_size, n_samples // 500)
    else:
        # For very large datasets, use larger minimum cluster size
        return max(base_min_size, n_samples // 1000)


def memory_efficient_hdbscan(embeddings, min_cluster_size=10):
    """Run HDBSCAN with strict memory optimization and better large dataset handling"""
    log_message("Running memory-efficient HDBSCAN...")
    n_samples = embeddings.shape[0]
    mem_before = log_memory()

    # Check available memory and adjust strategy accordingly
    available_memory = memory_limit_mb - (mem_before * 1.1)  # Leave 10% buffer

    if n_samples > 100000:
        log_message(f"Large dataset detected ({n_samples} samples)")
        log_message(f"Available memory: {available_memory:.0f} MB")

        # Try standard HDBSCAN first, but with strict memory monitoring
        try:
            log_message("Attempting memory-optimized HDBSCAN...")

            # Monitor memory during clustering
            if available_memory < 1000:  # Less than 1GB available
                log_message("Low memory available, using most conservative settings...")
                clusterer = hdbscan.HDBSCAN(
                    min_cluster_size=max(min_cluster_size, n_samples // 5000),  # Larger min size
                    approx_min_span_tree=True,
                    gen_min_span_tree=False,
                    core_dist_n_jobs=1,
                    cluster_selection_method='leaf'  # More memory efficient than 'eom'
                )
            else:
                clusterer = hdbscan.HDBSCAN(
                    min_cluster_size=min_cluster_size,
                    approx_min_span_tree=True,
                    gen_min_span_tree=False,
                    core_dist_n_jobs=1,
                    cluster_selection_method='eom'
                )

            # Monitor memory during fit_predict
            mem_during = log_memory()
            if mem_during > memory_limit_mb * 0.95:
                raise MemoryError("Memory usage too high during clustering")

            labels = clusterer.fit_predict(embeddings)

            # Validate results
            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            n_noise = list(labels).count(-1)

            log_message(f"Standard HDBSCAN completed: {n_clusters} clusters, {n_noise} noise points")
            log_memory()

            if n_clusters == 0:
                log_message("No clusters found, trying alternative approach...")
                raise ValueError("No clusters detected")

            return labels

        except (MemoryError, ValueError) as e:
            log_message(f"Standard HDBSCAN failed: {str(e)}")
            force_garbage_collection()

            # Fallback to memory-conservative sample-based approach
            return memory_conservative_clustering(embeddings, min_cluster_size, available_memory)

    else:
        # Standard HDBSCAN for smaller datasets
        log_message(f"Standard HDBSCAN for {n_samples} samples...")
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            cluster_selection_method='eom'
        )
        labels = clusterer.fit_predict(embeddings)

        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = list(labels).count(-1)
        log_message(f"HDBSCAN completed: {n_clusters} clusters, {n_noise} noise points")

        return labels


def memory_conservative_clustering(embeddings, min_cluster_size, available_memory_mb):
    """Ultra memory-conservative clustering for large datasets"""
    from sklearn.neighbors import NearestNeighbors
    from sklearn.cluster import DBSCAN

    n_samples = embeddings.shape[0]
    log_message(f"Starting memory-conservative clustering with {available_memory_mb:.0f} MB available")

    # Calculate sample size based on available memory
    # Rough estimate: HDBSCAN needs ~8 bytes per sample per sample for distance matrix
    # So for n samples, we need roughly n^2 * 8 bytes
    max_sample_for_memory = int(np.sqrt(available_memory_mb * 1024 * 1024 / 16))  # Conservative estimate
    sample_size = min(30000, max_sample_for_memory, n_samples // 4)  # Never more than 25% of data
    sample_size = max(sample_size, 5000)  # But at least 5000 for meaningful clustering

    log_message(f"Using conservative sample size: {sample_size}")

    # Use stratified sampling to get better representation
    sample_indices = np.random.choice(n_samples, sample_size, replace=False)
    sample_embeddings = embeddings[sample_indices]

    # Clear some memory
    force_garbage_collection()

    # Cluster the sample with conservative settings
    try:
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=max(min_cluster_size, sample_size // 1000),
            approx_min_span_tree=True,
            gen_min_span_tree=False,
            core_dist_n_jobs=1,
            cluster_selection_method='leaf'  # More memory efficient
        )
        sample_labels = clusterer.fit_predict(sample_embeddings)
        clustering_method = "HDBSCAN"

    except MemoryError:
        log_message("HDBSCAN on sample failed, using DBSCAN...")
        # Fallback to DBSCAN which is more memory efficient

        # Estimate eps efficiently
        nn = NearestNeighbors(n_neighbors=min(min_cluster_size, 20))  # Limit neighbors
        nn.fit(sample_embeddings[:5000])  # Use even smaller sample for eps estimation
        distances, _ = nn.kneighbors(sample_embeddings[:5000])
        eps = np.percentile(distances[:, -1], 85)

        # Clear NN from memory
        del nn, distances
        force_garbage_collection()

        dbscan = DBSCAN(eps=eps, min_samples=min_cluster_size)
        sample_labels = dbscan.fit_predict(sample_embeddings)
        clustering_method = "DBSCAN"

    n_sample_clusters = len(set(sample_labels)) - (1 if -1 in sample_labels else 0)
    log_message(f"Sample clustering ({clustering_method}) found {n_sample_clusters} clusters")

    if n_sample_clusters == 0:
        log_message("No clusters found in sample, returning all points as noise")
        return np.full(n_samples, -1, dtype=int)

    # Memory-efficient assignment of remaining points
    all_labels = np.full(n_samples, -1, dtype=int)
    all_labels[sample_indices] = sample_labels

    # Only assign labels to remaining points if we have clusters
    labeled_sample_mask = sample_labels != -1
    if np.any(labeled_sample_mask):
        labeled_points = sample_embeddings[labeled_sample_mask]
        labeled_labels = sample_labels[labeled_sample_mask]

        # Calculate cluster centroids (small memory footprint)
        unique_labels = np.unique(labeled_labels)
        cluster_centroids = np.zeros((len(unique_labels), embeddings.shape[1]))
        centroid_labels = []

        for i, label in enumerate(unique_labels):
            mask = labeled_labels == label
            cluster_centroids[i] = np.mean(labeled_points[mask], axis=0)
            centroid_labels.append(label)

        # Clear sample data from memory
        del sample_embeddings, labeled_points, labeled_labels
        force_garbage_collection()

        # Assign remaining points in very small batches to minimize memory usage
        unlabeled_mask = all_labels == -1
        unlabeled_indices = np.where(unlabeled_mask)[0]

        if len(unlabeled_indices) > 0:
            log_message(f"Assigning labels to {len(unlabeled_indices)} remaining points...")

            # Ultra-conservative batch size
            batch_size = min(5000, max(1000, available_memory_mb // 10))
            log_message(f"Using batch size: {batch_size}")

            for i in range(0, len(unlabeled_indices), batch_size):
                batch_indices = unlabeled_indices[i:i + batch_size]
                batch_points = embeddings[batch_indices]

                # Calculate distances to centroids (vectorized but memory-conscious)
                distances = np.linalg.norm(
                    batch_points[:, np.newaxis] - cluster_centroids[np.newaxis, :],
                    axis=2
                )

                # Conservative assignment - only assign closest points
                nearest_clusters = np.argmin(distances, axis=1)
                min_distances = np.min(distances, axis=1)

                # Use conservative threshold - only assign points very close to centroids
                threshold = np.percentile(min_distances, 50)  # Only closest 50%
                close_enough = min_distances < threshold

                batch_labels = np.full(len(batch_indices), -1)
                for j, (nearest_idx, close) in enumerate(zip(nearest_clusters, close_enough)):
                    if close:
                        batch_labels[j] = centroid_labels[nearest_idx]

                all_labels[batch_indices] = batch_labels

                # Clear batch data
                del batch_points, distances

                if (i // batch_size + 1) % 20 == 0:
                    log_message(f"Processed {i + batch_size}/{len(unlabeled_indices)} points")
                    if check_memory_limit():
                        force_garbage_collection()

    final_clusters = len(set(all_labels)) - (1 if -1 in all_labels else 0)
    final_noise = list(all_labels).count(-1)
    log_message(f"Memory-conservative clustering completed: {final_clusters} clusters, {final_noise} noise points")
    log_memory()

    return all_labels


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
    # First, try to load existing colors for current k
    existing_colors, _ = load_cluster_colors(barcode_dir, k)
    if existing_colors is not None:
        log_message(f"Using existing colors for k={k}")
        return existing_colors

    # Try to find any existing k-mer colors to use as base (k=6, k=7, or k=8)
    base_colors = None
    base_seq_clusters = None
    base_k = None

    for base_k_candidate in [6, 7, 8]:
        if base_k_candidate != k:  # Don't try to load from same k
            base_colors, base_seq_clusters = load_cluster_colors(barcode_dir, base_k_candidate)
            if base_colors is not None:
                base_k = base_k_candidate
                log_message(f"Found base colors from k={base_k} for k={k} color assignment")
                break

    # Generate colors for current clusters
    unique_clusters = sorted(c for c in set(labels) if c != -1)
    cmap = plt.get_cmap('tab20')
    colors = {}

    if base_colors is None:
        # No existing colors found, generate fresh colors
        log_message(f"Generating fresh colors for k={k}")
        colors = {cl: cmap(i % 20) for i, cl in enumerate(unique_clusters)}
        colors[-1] = (0.5, 0.5, 0.5, 0.5)  # grey for noise

        # Save these colors for future k-mer values
        save_cluster_colors(barcode_dir, k, colors, seq_ids, labels)
        return colors

    # Use base colors where sequences overlap
    current_seq_cluster_map = {seq_id: cluster for seq_id, cluster in zip(seq_ids, labels)}

    # Find which current clusters contain sequences from base clusters
    cluster_mapping = {}
    for seq_id, current_cluster in current_seq_cluster_map.items():
        if seq_id in base_seq_clusters:
            base_cluster = base_seq_clusters[seq_id]
            if current_cluster not in cluster_mapping:
                cluster_mapping[current_cluster] = {}
            if base_cluster not in cluster_mapping[current_cluster]:
                cluster_mapping[current_cluster][base_cluster] = 0
            cluster_mapping[current_cluster][base_cluster] += 1

    # Assign colors based on dominant base cluster in each current cluster
    used_base_colors = set()
    new_color_idx = len(base_colors)

    for current_cluster in unique_clusters:
        if current_cluster in cluster_mapping and cluster_mapping[current_cluster]:
            # Find the dominant base cluster
            dominant_base_cluster = max(cluster_mapping[current_cluster],
                                        key=cluster_mapping[current_cluster].get)
            if dominant_base_cluster in base_colors and dominant_base_cluster not in used_base_colors:
                colors[current_cluster] = base_colors[dominant_base_cluster]
                used_base_colors.add(dominant_base_cluster)
                log_message(
                    f"Cluster {current_cluster} inherits color from base cluster {dominant_base_cluster} (k={base_k})")
            else:
                # Assign new color
                colors[current_cluster] = cmap(new_color_idx % 20)
                new_color_idx += 1
                log_message(f"Cluster {current_cluster} gets new color (base color already used)")
        else:
            # New cluster with no base sequences, assign new color
            colors[current_cluster] = cmap(new_color_idx % 20)
            new_color_idx += 1
            log_message(f"Cluster {current_cluster} gets new color (no base overlap)")

    colors[-1] = (0.5, 0.5, 0.5, 0.5)  # grey for noise

    # Save colors for this k-mer value
    save_cluster_colors(barcode_dir, k, colors, seq_ids, labels)

    return colors


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

    n_sequences = vectors.shape[0]
    log_message(f"Loaded {n_sequences} vectors with {vectors.shape[1]} features")

    # UMAP
    if not os.path.exists(embeddings_file):
        embeddings = incremental_umap(vectors, embeddings_file)
    else:
        log_message(f"Found existing UMAP embeddings: {embeddings_file}")
        embeddings = np.load(embeddings_file)

    del vectors
    force_garbage_collection()

    # HDBSCAN with adaptive parameters
    if not os.path.exists(clusters_csv):
        # Use adaptive minimum cluster size
        adaptive_min_size = adaptive_min_cluster_size(n_sequences)
        log_message(f"Using adaptive min_cluster_size: {adaptive_min_size}")

        labels = memory_efficient_hdbscan(embeddings, adaptive_min_size)

        # Validate and log results
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = list(labels).count(-1)
        noise_pct = 100 * n_noise / len(labels)

        log_message(f"=== Clustering Results for {barcode_name} k={k} ===")
        log_message(f"Total sequences: {len(labels)}")
        log_message(f"Clusters found: {n_clusters}")
        log_message(f"Noise points: {n_noise} ({noise_pct:.1f}%)")

        if n_clusters > 0:
            cluster_sizes = [np.sum(labels == cl) for cl in set(labels) if cl != -1]
            log_message(
                f"Cluster size range: {min(cluster_sizes)}-{max(cluster_sizes)} (mean: {np.mean(cluster_sizes):.0f})")

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

    # Generate plot with memory-efficient colors
    colors = get_persistent_colors(barcode_dir, k, seq_ids, labels)

    # Create plot in memory-efficient way
    plt.figure(figsize=(12, 10))

    # Plot in batches to avoid memory issues with large datasets
    batch_size = 50000
    for i in range(0, len(embeddings), batch_size):
        end_idx = min(i + batch_size, len(embeddings))
        batch_embeddings = embeddings[i:end_idx]
        batch_labels = labels[i:end_idx]
        batch_colors = [colors.get(cl, (0.8, 0.8, 0.8, 0.5)) for cl in batch_labels]

        plt.scatter(batch_embeddings[:, 0], batch_embeddings[:, 1],
                    c=batch_colors, s=5, alpha=0.7)

    plt.title(f"UMAP + HDBSCAN ({barcode_name}, k={k})", fontsize=16)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()
    log_message(f"Saved plot: {plot_file}")

    force_garbage_collection()


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
    print("3 = Resume incomplete processing (check for existing files)")
    restart_mode = input("Enter 1, 2, or 3: ").strip()

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

    elif restart_mode == "3":
        # Smart resume: check what's already been processed
        log_message("Scanning for incomplete processing...")

        for barcode_idx, barcode_dir in enumerate(barcode_folders):
            barcode_name = os.path.basename(os.path.normpath(barcode_dir))

            for k in kmer_list:
                vector_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_vectors.npz")
                clusters_csv = os.path.join(barcode_dir, f"{barcode_name}_k{k}_clusters.csv")
                plot_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_UMAP_HDBSCAN.png")

                has_vectors = os.path.exists(vector_file)
                has_clusters = os.path.exists(clusters_csv)
                has_plot = os.path.exists(plot_file)

                status = []
                if has_vectors: status.append("vectors")
                if has_clusters: status.append("clusters")
                if has_plot: status.append("plot")

                if not (has_vectors and has_clusters and has_plot):
                    missing = []
                    if not has_vectors: missing.append("vectors")
                    if not has_clusters: missing.append("clusters")
                    if not has_plot: missing.append("plot")

                    print(f"{barcode_name} k={k}: INCOMPLETE - missing {', '.join(missing)}")
                else:
                    print(f"{barcode_name} k={k}: COMPLETE")

        confirm = input("\nProceed with smart resume? (y/n): ").strip().lower()
        if confirm != 'y':
            print("Resume cancelled.")
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

            # Check if this barcode/k combination is already complete (for smart resume)
            vector_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_vectors.npz")
            clusters_csv = os.path.join(barcode_dir, f"{barcode_name}_k{k}_clusters.csv")
            plot_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_UMAP_HDBSCAN.png")

            if restart_mode == "3" and all(os.path.exists(f) for f in [vector_file, clusters_csv, plot_file]):
                log_message(f"Skipping {barcode_name} k={k} (already complete)")
                continue

            log_message(f"Processing barcode {barcode_idx + 1}/{len(barcode_folders)}: {barcode_name} with k={k}")

            raw_dir = os.path.join(barcode_dir, "Raw")
            if not os.path.exists(raw_dir):
                log_message(f"No Raw folder in {barcode_name}, skipping.")
                continue

            fastq_files = glob.glob(os.path.join(raw_dir, "*_trimmed_cutadapt_*"))
            if not fastq_files:
                log_message(f"No FASTQ files in {barcode_name}, skipping.")
                continue

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
                log_message("Consider reducing memory_limit_mb or chunk sizes")
                continue
            except Exception as e:
                log_message(f"Unexpected error processing {barcode_name} k={k}: {str(e)}")
                continue

# Mode 2: Combined analysis (simplified version for memory efficiency)
elif mode == "2":
    log_message("Mode 2: Combined analysis not fully implemented in this memory-optimized version")
    log_message("Recommend using Mode 1 for memory-constrained systems")

else:
    log_message("Invalid mode selected")

end_time = datetime.datetime.now()
log_message(f"Processing finished. Total time: {end_time - start_time}")