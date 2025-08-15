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
import traceback

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
    try:
        process = psutil.Process(os.getpid())
        mem_mb = process.memory_info().rss / (1024 * 1024)
        log_message(f"RAM usage: {mem_mb:.2f} MB")
        return mem_mb
    except Exception as e:
        log_message(f"Error checking memory: {e}")
        return 0


def cleanup_memory():
    """Force garbage collection and memory cleanup"""
    gc.collect()
    log_message("Memory cleanup performed")


def kmer_vector(seq, k, kmer_index):
    seq = seq.upper()
    vec = np.zeros(len(kmer_index), dtype=int)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if all(b in bases for b in kmer):
            vec[kmer_index[kmer]] += 1
    return vec


def safe_vectorization(barcode_dir, barcode_name, k, fastq_files):
    """Safely perform vectorization with memory error handling"""
    try:
        vector_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_vectors.npz")

        log_message(f"Vectorizing {barcode_name} k={k}...")
        kmers = [''.join(p) for p in product(bases, repeat=k)]
        kmer_index = {km: idx for idx, km in enumerate(kmers)}

        seq_ids, file_labels, vectors = [], [], []

        for fpath in fastq_files:
            try:
                log_message(f"Reading {os.path.basename(fpath)}")
                file_vectors = []
                file_seq_ids = []
                file_labels_list = []

                for record in SeqIO.parse(fpath, "fastq"):
                    vec = kmer_vector(str(record.seq), k, kmer_index)
                    file_vectors.append(vec)
                    file_seq_ids.append(record.id)
                    file_labels_list.append(os.path.basename(fpath))

                # Add to main lists
                vectors.extend(file_vectors)
                seq_ids.extend(file_seq_ids)
                file_labels.extend(file_labels_list)

                # Clean up temporary lists
                del file_vectors, file_seq_ids, file_labels_list
                cleanup_memory()

                mem_usage = log_memory()
                # If memory usage is getting too high, force cleanup
                if mem_usage > 8000:  # 8GB threshold
                    cleanup_memory()

            except (MemoryError, Exception) as e:
                log_message(f"Error processing file {fpath}: {e}")
                log_message(f"Skipping file and continuing...")
                cleanup_memory()
                continue

        if not vectors:
            log_message(f"No vectors created for {barcode_name} k={k}")
            return False

        try:
            vectors = np.array(vectors, dtype=np.uint16)
            np.savez_compressed(vector_file, vectors=vectors, seq_ids=seq_ids, file_labels=file_labels)
            log_message(f"Saved {vectors.shape[0]} vectors to {vector_file}")

            # Clean up large arrays
            del vectors, seq_ids, file_labels
            cleanup_memory()
            return True

        except (MemoryError, Exception) as e:
            log_message(f"Error saving vectors for {barcode_name} k={k}: {e}")
            cleanup_memory()
            return False

    except (MemoryError, Exception) as e:
        log_message(f"Critical error in vectorization for {barcode_name} k={k}: {e}")
        log_message(f"Traceback: {traceback.format_exc()}")
        cleanup_memory()
        return False


def save_cluster_colors(barcode_dir, k, colors, seq_ids, labels):
    """Save cluster colors and sequence-to-cluster mapping to JSON file"""
    try:
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
        return True

    except Exception as e:
        log_message(f"Error saving cluster colors: {e}")
        return False


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
    try:
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

    except Exception as e:
        log_message(f"Error in get_persistent_colors: {e}")
        # Fallback to simple color assignment
        unique_clusters = sorted(c for c in set(labels) if c != -1)
        cmap = plt.get_cmap('tab20')
        colors = {cl: cmap(i % 20) for i, cl in enumerate(unique_clusters)}
        colors[-1] = (0.5, 0.5, 0.5, 0.5)  # grey for noise
        return colors


# ---------------------------
# Per-barcode analysis with error handling
# ---------------------------
def run_barcode_analysis(k, barcode_dir):
    """Run barcode analysis with comprehensive error handling"""
    barcode_name = os.path.basename(os.path.normpath(barcode_dir))

    try:
        vector_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_vectors.npz")
        clusters_csv = os.path.join(barcode_dir, f"{barcode_name}_k{k}_clusters.csv")
        embeddings_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_UMAP.npy")
        plot_file = os.path.join(barcode_dir, f"{barcode_name}_k{k}_UMAP_HDBSCAN.png")

        if not os.path.exists(vector_file):
            log_message(f"Vectors not found for {barcode_name} k={k}, skipping analysis.")
            return False

        log_message(f"Loading vectors for {barcode_name} k={k}...")
        try:
            data = np.load(vector_file, allow_pickle=True)
            vectors = data['vectors']
            seq_ids = data['seq_ids']
            file_labels = data['file_labels']
            log_memory()
        except (MemoryError, Exception) as e:
            log_message(f"Error loading vectors for {barcode_name} k={k}: {e}")
            return False

        # ---------------------------
        # UMAP
        # ---------------------------
        if not os.path.exists(embeddings_file):
            try:
                log_message(f"Running UMAP for {barcode_name} k={k}...")
                umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
                embeddings = umap_model.fit_transform(vectors)
                np.save(embeddings_file, embeddings)
                log_message(f"Saved UMAP embeddings: {embeddings_file}")
                cleanup_memory()
            except (MemoryError, Exception) as e:
                log_message(f"Error running UMAP for {barcode_name} k={k}: {e}")
                cleanup_memory()
                return False
        else:
            try:
                log_message(f"Found existing UMAP embeddings: {embeddings_file}")
                embeddings = np.load(embeddings_file)
            except (MemoryError, Exception) as e:
                log_message(f"Error loading UMAP embeddings for {barcode_name} k={k}: {e}")
                return False

        # ---------------------------
        # HDBSCAN
        # ---------------------------
        if not os.path.exists(clusters_csv):
            try:
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
                cleanup_memory()
            except (MemoryError, Exception) as e:
                log_message(f"Error running HDBSCAN for {barcode_name} k={k}: {e}")
                cleanup_memory()
                return False
        else:
            try:
                log_message(f"Found existing clusters CSV: {clusters_csv}")
                df_clusters = pd.read_csv(clusters_csv)
                labels = df_clusters['cluster'].values
            except Exception as e:
                log_message(f"Error loading clusters CSV for {barcode_name} k={k}: {e}")
                return False

        # ---------------------------
        # Get persistent colors and plot
        # ---------------------------
        try:
            colors = get_persistent_colors(barcode_dir, k, seq_ids, labels)
            point_colors = [colors.get(cl, (0.8, 0.8, 0.8, 0.5)) for cl in labels]

            # Plot
            plt.figure(figsize=(12, 10))
            plt.scatter(embeddings[:, 0], embeddings[:, 1], c=point_colors, s=5, alpha=0.7)
            plt.title(f"UMAP + HDBSCAN ({barcode_name}, k={k})", fontsize=16)
            plt.xlabel("UMAP1")
            plt.ylabel("UMAP2")
            plt.tight_layout()
            plt.savefig(plot_file, dpi=300)
            plt.close()
            log_message(f"Saved plot: {plot_file}")

            # Clean up
            del vectors, embeddings, labels, point_colors
            cleanup_memory()
            return True

        except Exception as e:
            log_message(f"Error creating plot for {barcode_name} k={k}: {e}")
            plt.close('all')  # Close any open plots
            cleanup_memory()
            return False

    except Exception as e:
        log_message(f"Unexpected error in barcode analysis for {barcode_name} k={k}: {e}")
        log_message(f"Traceback: {traceback.format_exc()}")
        plt.close('all')
        cleanup_memory()
        return False


# ---------------------------
# Combined analysis with error handling
# ---------------------------
def run_combined_analysis(k, barcode_dirs, base_colors=None):
    """Run combined analysis with comprehensive error handling"""
    try:
        log_message(f"Combined mode analysis for k={k}")
        combined_plot = os.path.join(log_dir, f"combined_k{k}_UMAP_HDBSCAN.png")
        combined_csv = os.path.join(log_dir, f"combined_k{k}_clusters.csv")
        embeddings_file = os.path.join(log_dir, f"combined_k{k}_embeddings.npy")
        temp_mmap_file = os.path.join(log_dir, f"combined_k{k}_vectors.dat")

        n_features = len([''.join(p) for p in product(bases, repeat=k)])
        total_vectors = 0

        # Count total vectors
        for barcode_dir in barcode_dirs:
            try:
                vector_file = glob.glob(os.path.join(barcode_dir, f"*k{k}_vectors.npz"))
                if vector_file:
                    data = np.load(vector_file[0], allow_pickle=True)
                    total_vectors += data['vectors'].shape[0]
                    del data  # Clean up immediately
            except Exception as e:
                log_message(f"Error counting vectors in {barcode_dir}: {e}")
                continue

        log_message(f"Total vectors found for k={k}: {total_vectors}")
        if total_vectors == 0:
            log_message("No vectors found. Skipping combined analysis.")
            return base_colors

        # ---------------------------
        # Load/Create memmap with error handling
        # ---------------------------
        try:
            if not os.path.exists(temp_mmap_file):
                log_message(f"Creating memmap file: {temp_mmap_file}")
                combined_vectors = np.memmap(temp_mmap_file, dtype='int32', mode='w+',
                                             shape=(total_vectors, n_features))
                all_seq_ids, all_file_labels = [], []
                idx = 0

                for barcode_dir in barcode_dirs:
                    try:
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
                            del data  # Clean up immediately
                            log_memory()
                    except Exception as e:
                        log_message(f"Error loading vectors from {barcode_dir}: {e}")
                        continue

                combined_vectors.flush()
            else:
                log_message(f"Using existing memmap: {temp_mmap_file}")
                combined_vectors = np.memmap(temp_mmap_file, dtype='int32', mode='r', shape=(total_vectors, n_features))
                all_seq_ids, all_file_labels = [], []

                for barcode_dir in barcode_dirs:
                    try:
                        vector_file = glob.glob(os.path.join(barcode_dir, f"*k{k}_vectors.npz"))
                        if vector_file:
                            data = np.load(vector_file[0], allow_pickle=True)
                            all_seq_ids.extend(data['seq_ids'])
                            all_file_labels.extend(data['file_labels'])
                            del data
                    except Exception as e:
                        log_message(f"Error loading sequence IDs from {barcode_dir}: {e}")
                        continue

        except (MemoryError, Exception) as e:
            log_message(f"Error creating/loading memmap: {e}")
            cleanup_memory()
            return base_colors

        # ---------------------------
        # UMAP with error handling
        # ---------------------------
        try:
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
                    try:
                        end = min(start + batch_size, total_vectors)
                        embeddings[start:end] = umap_model.transform(combined_vectors[start:end])
                        log_message(f"Transformed vectors {start} to {end}")
                        log_memory()
                    except (MemoryError, Exception) as e:
                        log_message(f"Error transforming batch {start}-{end}: {e}")
                        continue

                np.save(embeddings_file, embeddings)
                log_message(f"Saved embeddings: {embeddings_file}")
                cleanup_memory()
            else:
                log_message(f"Found existing embeddings: {embeddings_file}")
                embeddings = np.load(embeddings_file)

        except (MemoryError, Exception) as e:
            log_message(f"Error in UMAP processing: {e}")
            cleanup_memory()
            return base_colors

        # ---------------------------
        # HDBSCAN with error handling
        # ---------------------------
        try:
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
                cleanup_memory()
            else:
                log_message(f"Found existing CSV: {combined_csv}")
                df_clusters = pd.read_csv(combined_csv)
                labels = df_clusters['cluster'].values

        except (MemoryError, Exception) as e:
            log_message(f"Error in HDBSCAN processing: {e}")
            cleanup_memory()
            return base_colors

        # ---------------------------
        # Plotting with persistent colors
        # ---------------------------
        try:
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

            # Clean up
            del embeddings, labels, point_colors
            cleanup_memory()

            return base_colors

        except Exception as e:
            log_message(f"Error in plotting: {e}")
            plt.close('all')
            cleanup_memory()
            return base_colors

    except Exception as e:
        log_message(f"Unexpected error in combined analysis for k={k}: {e}")
        log_message(f"Traceback: {traceback.format_exc()}")
        plt.close('all')
        cleanup_memory()
        return base_colors


# ---------------------------
# Start with comprehensive error handling
# ---------------------------
log_message("Processing started.")
start_time = datetime.datetime.now()

print("Choose processing mode:")
print("1 = Process each barcode individually")
print("2 = Combine all barcodes into one analysis")
mode = input("Enter 1 or 2: ").strip()

barcode_folders = [d for d in glob.glob(os.path.join(wild2_dir, "*")) if os.path.isdir(d)]
log_message(f"Found {len(barcode_folders)} barcode folders")

# ---------------------------
# Mode 1: vectorization & per-barcode analysis
# ---------------------------
if mode == "1":
    kmer_list = [6, 7, 8]
    successful_barcodes = []
    failed_barcodes = []

    for k in kmer_list:
        log_message(f"Starting k={k} processing...")

        for barcode_dir in barcode_folders:
            barcode_name = os.path.basename(os.path.normpath(barcode_dir))
            log_message(f"Processing {barcode_name} with k={k}")

            try:
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
                    vectorization_success = True
                else:
                    vectorization_success = safe_vectorization(barcode_dir, barcode_name, k, fastq_files)

                if vectorization_success:
                    # Run per-barcode analysis
                    analysis_success = run_barcode_analysis(k, barcode_dir)
                    if analysis_success:
                        if barcode_name not in successful_barcodes:
                            successful_barcodes.append(barcode_name)
                        log_message(f"Successfully completed {barcode_name} k={k}")
                    else:
                        if barcode_name not in failed_barcodes:
                            failed_barcodes.append(barcode_name)
                        log_message(f"Analysis failed for {barcode_name} k={k}, but continuing...")
                else:
                    if barcode_name not in failed_barcodes:
                        failed_barcodes.append(barcode_name)
                    log_message(f"Vectorization failed for {barcode_name} k={k}, skipping to next barcode...")

            except Exception as e:
                log_message(f"Unexpected error processing {barcode_name} k={k}: {e}")
                log_message(f"Traceback: {traceback.format_exc()}")
                if barcode_name not in failed_barcodes:
                    failed_barcodes.append(barcode_name)
                log_message(f"Skipping {barcode_name} and continuing with next barcode...")
                cleanup_memory()
                continue

    # Summary
    log_message(f"Processing complete. Successful barcodes: {successful_barcodes}")
    log_message(f"Failed barcodes: {failed_barcodes}")

# ---------------------------
# Mode 2: combined analysis
# ---------------------------
elif mode == "2":
    kmer_list = [6, 7, 8]
    base_colors = None

    for k in kmer_list:
        try:
            log_message(f"Starting combined analysis for k={k}")
            base_colors = run_combined_analysis(k, barcode_folders, base_colors)
            log_message(f"Completed combined analysis for k={k}")
        except Exception as e:
            log_message(f"Failed combined analysis for k={k}: {e}")
            log_message(f"Traceback: {traceback.format_exc()}")
            log_message(f"Continuing with next k value...")
            cleanup_memory()
            continue

end_time = datetime.datetime.now()
log_message(f"Processing finished. Total time: {end_time - start_time}")
cleanup_memory()