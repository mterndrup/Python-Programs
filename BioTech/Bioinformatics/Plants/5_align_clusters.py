import os
import tempfile
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import subprocess
from collections import Counter, defaultdict
import numpy as np
import random
import time
from itertools import combinations
import hashlib

# Path to your MAFFT executable
MAFFT_PATH = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics"
              r"\Tools\MAFFT\mafft-7.526-win64-signed (1)\mafft-win\mafft.bat")

# Configuration for handling massive datasets
CLUSTER_SAMPLE_SIZE = 1000  # Sample size for clustering step
TREE_SAMPLE_SIZE = 200  # Maximum for tree construction
MIN_CLUSTER_SIZE = 5  # Minimum sequences per cluster to consider
SIMILARITY_THRESHOLD = 0.95  # For initial clustering


def sequence_hash(seq_str):
    """Create a hash for sequence deduplication."""
    return hashlib.md5(seq_str.encode()).hexdigest()[:16]


def deduplicate_sequences(records):
    """Remove exact duplicates while keeping track of counts."""
    unique_seqs = {}
    seq_counts = defaultdict(int)

    print("Deduplicating sequences...")
    for record in records:
        seq_str = str(record.seq).upper()
        seq_hash = sequence_hash(seq_str)

        if seq_hash not in unique_seqs:
            unique_seqs[seq_hash] = record
        seq_counts[seq_hash] += 1

    print(f"Deduplicated {len(records)} -> {len(unique_seqs)} unique sequences")
    return list(unique_seqs.values()), dict(seq_counts)


def cluster_by_length_and_composition(records, max_clusters=50):
    """Cluster sequences by length and basic composition for diversity sampling."""
    print("Clustering sequences by length and composition...")

    clusters = defaultdict(list)

    for record in records:
        seq_str = str(record.seq).upper()
        length = len(seq_str)

        # Basic composition signature
        gc_content = (seq_str.count('G') + seq_str.count('C')) / length if length > 0 else 0
        at_content = (seq_str.count('A') + seq_str.count('T')) / length if length > 0 else 0

        # Create cluster key based on length bins and composition
        length_bin = length // 10 * 10  # Group by 10bp bins
        gc_bin = int(gc_content * 10) / 10  # Group by 0.1 GC content bins

        cluster_key = (length_bin, gc_bin)
        clusters[cluster_key].append(record)

    print(f"Created {len(clusters)} clusters")

    # Sort clusters by size (largest first) and limit
    sorted_clusters = sorted(clusters.items(), key=lambda x: len(x[1]), reverse=True)
    return dict(sorted_clusters[:max_clusters])


def smart_sampling(clusters, target_size=CLUSTER_SAMPLE_SIZE):
    """Intelligently sample from clusters to maintain diversity."""
    total_sequences = sum(len(seqs) for seqs in clusters.values())
    print(f"Smart sampling from {total_sequences} sequences across {len(clusters)} clusters")

    if total_sequences <= target_size:
        # Use all sequences if under target
        all_seqs = []
        for seqs in clusters.values():
            all_seqs.extend(seqs)
        return all_seqs

    sampled_sequences = []
    remaining_budget = target_size
    min_per_cluster = max(1, target_size // (len(clusters) * 4))

    cluster_samples = {}
    used_seq_ids = set()

    for cluster_key, sequences in clusters.items():
        if len(sequences) >= MIN_CLUSTER_SIZE and remaining_budget > 0:
            proportion = len(sequences) / total_sequences
            target_from_cluster = max(min_per_cluster, int(proportion * target_size))
            actual_sample = min(target_from_cluster, len(sequences), remaining_budget)

            if actual_sample > 0:
                sampled = random.sample(sequences, actual_sample)
                cluster_samples[cluster_key] = sampled
                for seq in sampled:
                    used_seq_ids.add(seq.id)
                remaining_budget -= actual_sample

    if remaining_budget > 0:
        large_clusters = sorted([(k, v) for k, v in clusters.items()
                                 if len(v) >= MIN_CLUSTER_SIZE and k in cluster_samples],
                                key=lambda x: len(x[1]), reverse=True)

        for cluster_key, sequences in large_clusters:
            if remaining_budget <= 0:
                break

            current_sample_size = len(cluster_samples.get(cluster_key, []))
            additional_sample = min(remaining_budget,
                                    len(sequences) - current_sample_size,
                                    len(sequences) // 4)

            if additional_sample > 0:
                available_seqs = [s for s in sequences if s.id not in used_seq_ids]
                if available_seqs and additional_sample > 0:
                    additional_sample = min(additional_sample, len(available_seqs))
                    additional_seqs = random.sample(available_seqs, additional_sample)
                    cluster_samples[cluster_key].extend(additional_seqs)
                    for seq in additional_seqs:
                        used_seq_ids.add(seq.id)
                    remaining_budget -= additional_sample

    for samples in cluster_samples.values():
        sampled_sequences.extend(samples)

    print(f"Sampled {len(sampled_sequences)} sequences for analysis")
    return sampled_sequences


def filter_top_length_sequences(records):
    """Filter sequences to top-length (most common length)."""
    lengths = [len(rec.seq) for rec in records]
    length_counts = Counter(lengths)
    top_length = length_counts.most_common(1)[0][0]
    filtered = [rec for rec in records if len(rec.seq) == top_length]
    return filtered, top_length, len(records)


def run_mafft_progressive(records):
    """Run MAFFT with progressive alignment for large datasets."""
    print(f"Progressive alignment of {len(records)} sequences...")
    start_time = time.time()

    with tempfile.TemporaryDirectory() as tmpdir:
        input_fasta = os.path.join(tmpdir, "input.fasta")
        output_fasta = os.path.join(tmpdir, "aligned.fasta")
        SeqIO.write(records, input_fasta, "fasta")

        try:
            if len(records) > 500:
                cmd = [MAFFT_PATH, "--retree", "1", "--maxiterate", "0", "--quiet", input_fasta]
                timeout = 600
            elif len(records) > 100:
                cmd = [MAFFT_PATH, "--retree", "2", "--maxiterate", "1", "--quiet", input_fasta]
                timeout = 300
            else:
                cmd = [MAFFT_PATH, "--auto", "--quiet", input_fasta]
                timeout = 180

            print(f"Running MAFFT: {' '.join(cmd[1:-2])}")

            subprocess.run(cmd, stdout=open(output_fasta, "w"), check=True, timeout=timeout)

        except subprocess.TimeoutExpired:
            print(f"MAFFT alignment timed out after {timeout / 60:.1f} minutes")
            return None
        except FileNotFoundError:
            print("MAFFT not found. Please check MAFFT_PATH.")
            return None
        except subprocess.CalledProcessError as e:
            print(f"MAFFT failed with error: {e}")
            return None

        aligned = list(SeqIO.parse(output_fasta, "fasta"))
        elapsed = time.time() - start_time
        print(f"Alignment completed in {elapsed:.1f} seconds")
        return aligned


def build_tree_subset(alignment, max_seqs_for_tree=TREE_SAMPLE_SIZE, output_file=None):
    """Build tree from a subset if alignment is too large."""
    print(f"Building phylogenetic tree...")

    alignment_records = list(alignment)

    if len(alignment_records) > max_seqs_for_tree:
        print(f"Subsampling {len(alignment_records)} -> {max_seqs_for_tree} sequences for tree construction")
        indices = np.linspace(0, len(alignment_records) - 1, max_seqs_for_tree, dtype=int)
        indices = [int(i) for i in indices]
        tree_records = [alignment_records[i] for i in indices]
        tree_alignment = MultipleSeqAlignment(tree_records)
    else:
        tree_records = alignment_records
        tree_alignment = alignment

    start_time = time.time()

    try:
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(tree_alignment)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

        if output_file:
            Phylo.write(tree, output_file, "newick")
            print(f"Tree written to {output_file}")

        elapsed = time.time() - start_time
        print(f"Tree construction completed in {elapsed:.1f} seconds")
        return tree, dm, tree_records

    except Exception as e:
        print(f"Tree construction failed: {e}")
        return None, None, None


def analyze_large_dataset_centrality(dm, tree_records, all_aligned_records, original_records,
                                     log_file, blast_fasta, fasta_path, total_original, region):
    """Analyze centrality for large datasets - select only central/representative sequences."""
    ids = dm.names
    dist_matrix = np.array(dm)
    avg_distances = dist_matrix.mean(axis=1)

    n_central = min(20, len(ids) // 2, len(ids))
    central_indices = np.argsort(avg_distances)[:n_central]
    divergent_indices = np.argsort(avg_distances)[-10:]

    tree_record_map = {rec.id: rec for rec in tree_records}
    original_map = {rec.id: rec for rec in original_records}

    sample_size = min(50, len(all_aligned_records))
    sample_records = random.sample(list(all_aligned_records), sample_size)

    consensus_seq = ""
    if sample_records:
        for i in range(len(sample_records[0].seq)):
            column = [rec.seq[i] for rec in sample_records]
            consensus_seq += Counter(column).most_common(1)[0][0]

        differences = [sum(1 for a, b in zip(rec.seq, consensus_seq) if a != b) for rec in sample_records]
    else:
        differences = [0]

    with open(log_file, "w") as f:
        f.write(f"{region} Cluster QC Report - Large Dataset Analysis\n")
        f.write(f"FASTA: {fasta_path}\n")
        f.write(f"Total sequences in file: {total_original}\n")
        f.write(f"Sequences used in alignment: {len(all_aligned_records)}\n")
        f.write(f"Sequences used in tree: {len(tree_records)}\n")
        f.write(f"Alignment length: {len(all_aligned_records[0].seq) if all_aligned_records else 'N/A'} bp\n")
        f.write(f"Consensus calculated from: {len(sample_records)} sequences\n")
        if differences:
            f.write(f"Average differences from consensus: {np.mean(differences):.2f}\n")
            f.write(f"Max differences: {np.max(differences)}\n")
            f.write(f"Min differences: {np.min(differences)}\n")
        f.write(f"\nNote: Analysis performed on representative subset due to large dataset size\n")
        f.write(f"Selection strategy: Central sequences only (divergent sequences excluded as likely noise)\n\n")

        f.write(f"{n_central} Most Central/Representative Sequences (SELECTED FOR BLAST):\n")
        for i in central_indices:
            f.write(f"{ids[i]};avg_dist={avg_distances[i]:.4f}\n")

        f.write(f"\n10 Most Divergent Sequences (NOT SELECTED - likely noise/artifacts):\n")
        for i in divergent_indices:
            f.write(f"{ids[i]};avg_dist={avg_distances[i]:.4f}\n")

    selected_records = []
    selected_ids = set()

    for i in central_indices:
        seq_id = ids[i]
        if seq_id in original_map and seq_id not in selected_ids:
            selected_records.append(original_map[seq_id])
            selected_ids.add(seq_id)

    SeqIO.write(selected_records, blast_fasta, "fasta")

    print(f"QC log written to {log_file}")
    print(f"BLAST-ready FASTA written to {blast_fasta} ({len(selected_records)} central sequences only)")
    print(f"Excluded {len(divergent_indices)} highly divergent sequences (likely noise)")


def process_sample(barcode_num, base_path, region):
    if region == "ITS2":
        fasta_path = os.path.join(base_path, f"barcode{barcode_num}",
                                  "ITSx-All", f"barcode{barcode_num}_combined_ITSx.ITS2_centroids.fasta")
    elif region == "LSU":
        fasta_path = os.path.join(base_path, f"barcode{barcode_num}",
                                  "ITSx-All", f"barcode{barcode_num}_combined_ITSx.LSU_centroids.fasta")
    else:
        print(f"Unknown region: {region}")
        return

    if not os.path.isfile(fasta_path):
        print(f"FASTA file not found:\n{fasta_path}")
        return

    print(f"Loading {region} sequences from {fasta_path}")
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        print("No sequences found in FASTA!")
        return

    total_original = len(records)
    print(f"Loaded {total_original:,} {region} sequences")

    if total_original < 10:
        print("Too few sequences for meaningful analysis")
        return

    unique_records, seq_counts = deduplicate_sequences(records)
    filtered_records, top_length, _ = filter_top_length_sequences(unique_records)
    print(f"After length filtering: {len(filtered_records):,} sequences of length {top_length}")

    if len(filtered_records) < 10:
        print("Too few sequences after filtering")
        return

    if len(filtered_records) > CLUSTER_SAMPLE_SIZE:
        print(f"Large dataset detected ({len(filtered_records):,} sequences)")
        clusters = cluster_by_length_and_composition(filtered_records)
        sampled_records = smart_sampling(clusters, CLUSTER_SAMPLE_SIZE)
    else:
        sampled_records = filtered_records

    print(f"Using {len(sampled_records)} {region} sequences for alignment")

    aligned_records = run_mafft_progressive(sampled_records)
    if aligned_records is None:
        print("MAFFT alignment failed!")
        return

    log_dir = os.path.join(base_path, f"barcode{barcode_num}", "Logs")
    os.makedirs(log_dir, exist_ok=True)

    tree_file = os.path.join(log_dir, f"barcode{barcode_num}_{region}_clusters_tree.nwk")
    blast_fasta = os.path.join(log_dir, f"barcode{barcode_num}_{region}_clusters_selected.fasta")
    log_file = os.path.join(log_dir, f"barcode{barcode_num}_{region}_clusters_QC.log")

    tree, dm, tree_records = build_tree_subset(MultipleSeqAlignment(aligned_records),
                                               output_file=tree_file)

    if tree is None or dm is None:
        print("Tree construction failed, skipping centrality analysis")
        return None

    analyze_large_dataset_centrality(dm, tree_records, aligned_records, filtered_records,
                                     log_file, blast_fasta, fasta_path, total_original, region)

    return tree_file, log_file, blast_fasta


def main():
    base_path = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus-All"

    print("Large-Scale Bioinformatics Pipeline")
    print(f"Configuration for massive datasets:")
    print(f"- Clustering sample size: {CLUSTER_SAMPLE_SIZE:,}")
    print(f"- Tree construction limit: {TREE_SAMPLE_SIZE:,}")
    print(f"- Minimum cluster size: {MIN_CLUSTER_SIZE}")
    print()

    region = input("Process which region? (ITS2/LSU): ").strip().upper()
    if region not in ["ITS2", "LSU"]:
        print("Invalid region, please enter ITS2 or LSU.")
        return

    choice = input("Process one barcode or all? (one/all): ").strip().lower()
    if choice == "one":
        barcode_num = input("Enter barcode number (e.g., 73): ").strip()
        if not barcode_num.isdigit():
            print("Invalid barcode number")
            return
        barcode_num = barcode_num.zfill(2)
        print(f"Processing barcode {barcode_num} ({region})...")
        start_time = time.time()
        result = process_sample(barcode_num, base_path, region)
        elapsed = time.time() - start_time

        if result:
            tree_file, log_file, blast_fasta = result
            print(f"\nCompleted in {elapsed / 60:.1f} minutes")
            print(f"Tree saved to: {tree_file}")
            print(f"QC log saved to: {log_file}")
            print(f"BLAST FASTA saved to: {blast_fasta}")
        else:
            print("Error processing sample.")

    elif choice == "all":
        total_start = time.time()
        successful = 0
        failed = 0

        for name in os.listdir(base_path):
            if name.startswith("barcode") and os.path.isdir(os.path.join(base_path, name)):
                barcode_num = name.replace("barcode", "")
                print(f"\n{'=' * 50}")
                print(f"Processing barcode {barcode_num} ({region})")
                print(f"{'=' * 50}")
                start_time = time.time()
                result = process_sample(barcode_num, base_path, region)
                elapsed = time.time() - start_time

                if result:
                    print(f"✓ Completed barcode {barcode_num} in {elapsed / 60:.1f} minutes")
                    successful += 1
                else:
                    print(f"✗ Error processing barcode {barcode_num}")
                    failed += 1

        total_elapsed = time.time() - total_start
        print(f"\n{'=' * 50}")
        print(f"FINAL SUMMARY ({region})")
        print(f"{'=' * 50}")
        print(f"Successful: {successful}")
        print(f"Failed: {failed}")
        print(f"Total time: {total_elapsed / 3600:.1f} hours")

    else:
        print("Invalid choice, please enter 'one' or 'all'.")


if __name__ == "__main__":
    main()
