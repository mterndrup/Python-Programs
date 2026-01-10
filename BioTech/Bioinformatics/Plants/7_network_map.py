import os
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# === User configuration ===
BASE_PATH = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus-All"
FASTA_FILE = os.path.join(BASE_PATH, "Cactus_All_top3_ITS_LSU_combined.fasta")
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
        # Assumes sample ID is everything before first underscore
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
print(f"Alignment complete: {len(alignment)} sequences, length {alignment[0].seq.__len__()}")

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
G = nx.Graph()
for rec in alignment:
    G.add_node(rec.id)

# Connect sequences that are "close" (distance threshold)
threshold = 0.05  # adjust as needed
n = len(alignment)
for i in range(n):
    for j in range(i+1, n):
        if dm[i,j] <= threshold:
            G.add_edge(alignment[i].id, alignment[j].id, weight=dm[i,j])

# === Step 7: Plot network ===
plt.figure(figsize=(12,12))
pos = nx.spring_layout(G, k=0.5, seed=42)
nx.draw_networkx_nodes(G, pos, node_color='skyblue', node_size=500)
nx.draw_networkx_edges(G, pos, alpha=0.5)
nx.draw_networkx_labels(G, pos, font_size=8)
plt.title("Cactus Lineage Network (ITS2+LSU)")
plt.axis('off')
plt.tight_layout()
plt.show()

# Optional: save figure
#plt.savefig(os.path.join(BASE_PATH, "Cactus_Lineage_Network.png"), dpi=300)
#print("Network saved as Cactus_Lineage_Network.png")
