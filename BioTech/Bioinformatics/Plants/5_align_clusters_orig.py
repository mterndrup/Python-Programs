import os
import tempfile
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import subprocess
from collections import Counter
import numpy as np

# Path to your MAFFT executable
MAFFT_PATH = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics"
              r"\Tools\MAFFT\mafft-7.526-win64-signed (1)\mafft-win\mafft.bat")

def filter_top_length_sequences(records):
    """Filter sequences to top-length (most common length)."""
    lengths = [len(rec.seq) for rec in records]
    length_counts = Counter(lengths)
    top_length = length_counts.most_common(1)[0][0]
    filtered = [rec for rec in records if len(rec.seq) == top_length]
    return filtered, top_length, len(records)

def run_mafft(records):
    """Run MAFFT alignment and return aligned SeqRecords."""
    with tempfile.TemporaryDirectory() as tmpdir:
        input_fasta = os.path.join(tmpdir, "input.fasta")
        output_fasta = os.path.join(tmpdir, "aligned.fasta")
        SeqIO.write(records, input_fasta, "fasta")
        try:
            subprocess.run([MAFFT_PATH, "--auto", "--quiet", input_fasta],
                           stdout=open(output_fasta, "w"), check=True)
        except FileNotFoundError:
            print("MAFFT not found. Please check MAFFT_PATH.")
            return None
        aligned = list(SeqIO.parse(output_fasta, "fasta"))
        return aligned

def build_tree(alignment, output_file=None):
    """Build a neighbor-joining tree from aligned sequences."""
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    if output_file:
        Phylo.write(tree, output_file, "newick")
        print(f"Tree written to {output_file}")
    return tree, dm

def analyze_centrality(dm, aligned_records, original_records, log_file, blast_fasta, fasta_path):
    """Compute central/divergent sequences, log QC, and write BLAST FASTA."""
    ids = dm.names
    dist_matrix = np.array(dm)
    avg_distances = dist_matrix.mean(axis=1)

    central_indices = np.argsort(avg_distances)[:5]
    divergent_indices = np.argsort(avg_distances)[-10:]

    original_map = {rec.id: rec for rec in original_records}

    # Compute consensus info
    consensus_seq = ""
    for i in range(len(aligned_records[0].seq)):
        column = [rec.seq[i] for rec in aligned_records]
        consensus_seq += Counter(column).most_common(1)[0][0]
    differences = [sum(1 for a,b in zip(rec.seq, consensus_seq) if a != b) for rec in aligned_records]

    # Write QC log
    with open(log_file, "w") as f:
        f.write("ITS2 Cluster QC Report\n")
        f.write(f"FASTA: {fasta_path}\n")
        f.write(f"Total sequences: {len(original_records)}\n")
        f.write(f"Sequences of top length ({len(aligned_records[0].seq)} bp): {len(aligned_records)}\n")
        f.write(f"Consensus length: {len(consensus_seq)}\n")
        f.write(f"Average differences from consensus: {np.mean(differences):.2f}\n")
        f.write(f"Max differences: {np.max(differences)}\n")
        f.write(f"Min differences: {np.min(differences)}\n\n")

        f.write("5 Most Central Sequences:\n")
        for i in central_indices:
            f.write(f"{ids[i]};avg_dist={avg_distances[i]:.4f}\n")

        f.write("\n10 Most Divergent Sequences:\n")
        for i in divergent_indices:
            f.write(f"{ids[i]};avg_dist={avg_distances[i]:.4f}\n")

    # Write BLAST-ready FASTA
    selected_records = []
    for i in central_indices:
        sid = ids[i]
        if sid in original_map:
            selected_records.append(original_map[sid])
    for i in divergent_indices:
        sid = ids[i]
        if sid in original_map:
            selected_records.append(original_map[sid])
    SeqIO.write(selected_records, blast_fasta, "fasta")

    print(f"QC log written to {log_file}")
    print(f"BLAST-ready FASTA written to {blast_fasta}")

def process_sample(barcode_num, base_path):
    fasta_path = os.path.join(base_path, f"barcode{barcode_num}",
                              "ITSx-All", f"barcode{barcode_num}_combined_ITSx.ITS2_centroids.fasta")
    if not os.path.isfile(fasta_path):
        print(f"FASTA file not found:\n{fasta_path}")
        return

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        print("No sequences found in FASTA!")
        return

    filtered_records, top_length, total_sequences = filter_top_length_sequences(records)
    aligned_records = run_mafft(filtered_records)
    if aligned_records is None:
        print("MAFFT alignment failed!")
        return

    log_dir = os.path.join(base_path, f"barcode{barcode_num}", "Logs")
    os.makedirs(log_dir, exist_ok=True)

    tree_file = os.path.join(log_dir, f"barcode{barcode_num}_ITS2_clusters_tree.nwk")
    blast_fasta = os.path.join(log_dir, f"barcode{barcode_num}_ITS2_clusters_selected.fasta")
    log_file = os.path.join(log_dir, f"barcode{barcode_num}_ITS2_clusters_QC.log")

    tree, dm = build_tree(MultipleSeqAlignment(aligned_records), output_file=tree_file)
    analyze_centrality(dm, aligned_records, filtered_records, log_file, blast_fasta, fasta_path)

    return tree_file, log_file, blast_fasta

def main():
    base_path = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus2"

    choice = input("Process one barcode or all? (one/all): ").strip().lower()
    if choice == "one":
        barcode_num = input("Enter barcode number (e.g., 73): ").strip()
        if not barcode_num.isdigit():
            print("Invalid barcode number")
            return
        barcode_num = barcode_num.zfill(2)
        print(f"Processing barcode {barcode_num}...")
        result = process_sample(barcode_num, base_path)
        if result:
            tree_file, log_file, blast_fasta = result
            print(f"Tree saved to: {tree_file}")
            print(f"QC log saved to: {log_file}")
            print(f"BLAST FASTA saved to: {blast_fasta}")
        else:
            print("Error processing sample.")
    elif choice == "all":
        for name in os.listdir(base_path):
            if name.startswith("barcode") and os.path.isdir(os.path.join(base_path, name)):
                barcode_num = name.replace("barcode", "")
                print(f"Processing barcode {barcode_num}...")
                result = process_sample(barcode_num, base_path)
                if result:
                    tree_file, log_file, blast_fasta = result
                    print(f"Tree saved to: {tree_file}")
                    print(f"QC log saved to: {log_file}")
                    print(f"BLAST FASTA saved to: {blast_fasta}")
                else:
                    print(f"Error processing barcode {barcode_num}.")
    else:
        print("Invalid choice, please enter 'one' or 'all'.")

if __name__ == "__main__":
    main()
