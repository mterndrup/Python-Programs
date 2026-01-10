import os
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# === User configuration ===
BASE_PATH = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus-All"
FASTA_FILE = os.path.join(BASE_PATH, "Cactus_All_top3_ITS_LSU_combined_select.fasta")
MAFFT_PATH = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Tools\MAFFT\mafft-7.526-win64-signed (1)\mafft-win\mafft.bat"

# Ask user whether to keep all sequences or one per sample
choice = input("Show all nodes or one per sample? (all/one): ").strip().lower()
keep_all = choice == "all"

# === Step 1: Read sequences ===
records = list(SeqIO.parse(FASTA_FILE, "fasta"))
print(f"Loaded {len(records)} sequences from {FASTA_FILE}")

# Optionally reduce to one sequence per sample
if not keep_all:
    unique_samples = {}
    for rec in records:
        sample_id = rec.id.split("_")[0]
        if sample_id not in unique_samples:
            unique_samples[sample_id] = rec
    records = list(unique_samples.values())
    print(f"Reduced to {len(records)} sequences (one per sample)")

# === Step 2: Write sequences to temp FASTA for MAFFT ===
aligned_input = os.path.join(BASE_PATH, "temp_input.fasta")
SeqIO.write(records, aligned_input, "fasta")

# === Step 3: Run MAFFT ===
aligned_output = os.path.join(BASE_PATH, "aligned_sequences.fasta")
print("Running MAFFT alignment...")
mafft_cline = MafftCommandline(MAFFT_PATH, input=aligned_input, quiet=True)
stdout, stderr = mafft_cline()
with open(aligned_output, "w") as handle:
    handle.write(stdout)

# === Step 4: Read aligned sequences ===
alignment = list(AlignIO.read(aligned_output, "fasta"))
print(f"Alignment complete: {len(alignment)} sequences, length {len(alignment[0].seq)}")

# === Step 5: Build distance matrix ===
def build_distance_matrix(aligned_records):
    n = len(aligned_records)
    dm = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            seq1 = str(aligned_records[i].seq)
            seq2 = str(aligned_records[j].seq)
            matches = sum(a==b for a,b in zip(seq1, seq2))
            distance = 1 - matches / len(seq1)
            dm[i,j] = dm[j,i] = distance
    return dm

dm = build_distance_matrix(alignment)

# === Step 6: Build network graph ===
G = nx.DiGraph()  # use directed graph for hybrid arrows
for rec in alignment:
    G.add_node(rec.id)

# Connect sequences below threshold
threshold = 0.05  # adjust as needed
n = len(alignment)
for i in range(n):
    for j in range(i+1, n):
        if dm[i,j] <= threshold:
            # undirected connection for similar sequences
            G.add_edge(alignment[i].id, alignment[j].id)
            G.add_edge(alignment[j].id, alignment[i].id)

# === Step 7: Add hybrid nodes (optional) ===
# For demonstration, we'll create hybrids if a node is "between" two parents
# Here, a simple heuristic: if a sequence is within threshold of two distinct parents, consider hybrid
hybrids = []
for i in range(n):
    close_parents = [j for j in range(n) if 0 < dm[i,j] <= threshold]
    if len(close_parents) >= 2:
        parents = [alignment[j].id for j in close_parents[:2]]
        hybrid_name = f"Hybrid_{alignment[i].id}"
        G.add_node(hybrid_name)
        # arrows from parents to hybrid
        for p in parents:
            G.add_edge(p, hybrid_name)
        hybrids.append(hybrid_name)

# === Step 8: Plot network ===
plt.figure(figsize=(12,12))
pos = nx.spring_layout(G, k=0.5, seed=42)

# separate hybrid nodes for coloring
node_colors = ['lightgreen' if n in hybrids else 'skyblue' for n in G.nodes]

nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500)
nx.draw_networkx_edges(G, pos, alpha=0.5, arrowsize=15)
nx.draw_networkx_labels(G, pos, font_size=8)
plt.title("Cactus Lineage Network with Hybrids (ITS2+LSU)")
plt.axis('off')
plt.tight_layout()
plt.show()

plt.savefig(os.path.join(BASE_PATH, "Cactus_Lineage_Network_Hybrids.png"), dpi=300)
print("Network saved as Cactus_Lineage_Network_Hybrids.png")
