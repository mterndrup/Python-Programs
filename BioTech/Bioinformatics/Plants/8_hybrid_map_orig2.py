import os
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.manifold import MDS  # for similarity-based layout

# === User configuration ===
BASE_PATH = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus-All"
FASTA_FILE = os.path.join(BASE_PATH, "Cactus_All_top3_ITS_LSU_combined_select.fasta")
MAFFT_PATH = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Tools\MAFFT\mafft-7.526-win64-signed (1)\mafft-win\mafft.bat"

# === Step 0: User options ===
mode = input("Build network by threshold or MSN? (threshold/msn): ").strip().lower()
choice = input("Show all nodes or one per sample? (all/one): ").strip().lower()
keep_all = choice == "all"

if mode == "threshold":
    layout_choice = input("Threshold layout: similarity-based or evenly spaced? (similarity/even): ").strip().lower()
else:
    layout_choice = "similarity"  # MSN uses similarity-based layout

# === Step 1: Read sequences ===
records = list(SeqIO.parse(FASTA_FILE, "fasta"))
print(f"Loaded {len(records)} sequences from {FASTA_FILE}")

if not keep_all:
    unique_samples = {}
    for rec in records:
        sample_id = rec.id.split("_")[0]
        if sample_id not in unique_samples:
            unique_samples[sample_id] = rec
    records = list(unique_samples.values())
    print(f"Reduced to {len(records)} sequences (one per sample)")
else:
    print("Keeping all sequences per sample.")

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
        for j in range(i + 1, n):
            seq1 = str(aligned_records[i].seq)
            seq2 = str(aligned_records[j].seq)
            matches = sum(a == b for a, b in zip(seq1, seq2))
            distance = 1 - matches / len(seq1)
            dm[i, j] = dm[j, i] = distance
    return dm

dm = build_distance_matrix(alignment)

# === Step 6: Build network graph ===
threshold = 0.05
hybrid_threshold = 0.6  # 60% likelihood threshold for calling something a hybrid
n = len(alignment)
G = nx.Graph() if mode == "msn" else nx.DiGraph()  # MSN = undirected, threshold = directed

for rec in alignment:
    G.add_node(rec.id)

if mode == "threshold":
    for i in range(n):
        for j in range(i + 1, n):
            if dm[i, j] <= threshold:
                G.add_edge(alignment[i].id, alignment[j].id, weight=1 - dm[i, j])
                G.add_edge(alignment[j].id, alignment[i].id, weight=1 - dm[i, j])
else:  # MSN
    G_full = nx.Graph()
    for i in range(n):
        for j in range(i + 1, n):
            G_full.add_edge(alignment[i].id, alignment[j].id, weight=dm[i,j])
    G = nx.minimum_spanning_tree(G_full)
    print(f"MSN built with {len(G.nodes)} nodes and {len(G.edges)} edges")

# === Step 7: Hybrid likelihood calculation ===
def calculate_hybrid_likelihood(seq_idx, distance_matrix, threshold=0.05):
    distances = distance_matrix[seq_idx]
    n_seqs = len(distances)
    potential_parents = [(j, distances[j]) for j in range(n_seqs) if j != seq_idx and 0 < distances[j] <= threshold]
    if len(potential_parents) < 2:
        return 0.0, []
    potential_parents.sort(key=lambda x: x[1])
    factors = []
    n_close = len(potential_parents)
    factors.append(min(n_close / 3.0, 1.0))
    if n_close >= 2:
        dist_to_parents = [d for _, d in potential_parents[:3]]
        avg_dist = np.mean(dist_to_parents)
        var_dist = np.var(dist_to_parents)
        factors.extend([max(0, 1 - abs(avg_dist - threshold/2)/(threshold/2)), 1 - min(var_dist/(threshold/4),1.0)])
    if len(potential_parents) >= 2:
        closest_dist = potential_parents[0][1]
        second_dist = potential_parents[1][1]
        factors.append(closest_dist / second_dist if second_dist > 0 else 0)
    weights = [0.3,0.25,0.25,0.2][:len(factors)]
    weights = [w/sum(weights) for w in weights]
    likelihood = sum(f*w for f,w in zip(factors, weights))
    likelihood = min(likelihood, 1.0)
    top_parents = potential_parents[:2]
    return likelihood, top_parents

# Compute hybrid data
hybrid_data = {}
hybrids = []
nodes_to_remove = []

for i in range(n):
    likelihood, parents = calculate_hybrid_likelihood(i, dm, threshold)
    seq_id = alignment[i].id
    if likelihood > 0:
        hybrid_data[seq_id] = {'likelihood': likelihood, 'parents': parents, 'percentage': likelihood*100}
    if likelihood >= hybrid_threshold and len(parents) >= 2:
        sample_name = seq_id.split('_')[0]
        hybrid_name = f"Hybrid_{sample_name}"
        G.add_node(hybrid_name)
        for parent_idx, distance in parents:
            G.add_edge(alignment[parent_idx].id, hybrid_name, weight=1 - distance)
        hybrids.append(hybrid_name)
        hybrid_data[hybrid_name] = hybrid_data[seq_id]
        nodes_to_remove.append(seq_id)

for node in nodes_to_remove:
    if G.has_node(node):
        G.remove_node(node)
        print(f"Removed original node {node}, replaced with hybrid version")

# === Helper functions for labels/colors ===
def clean_label(node_id):
    if node_id.startswith("Hybrid_"):
        base_name = node_id.replace("Hybrid_", "")
        if node_id in hybrid_data:
            percentage = hybrid_data[node_id]['percentage']
            return f"{base_name}\n({percentage:.1f}%)"
        return base_name
    parts = node_id.split("_")
    if len(parts) >= 2:
        sample = parts[0]
        if node_id in hybrid_data and hybrid_data[node_id]['percentage']>30:
            percentage = hybrid_data[node_id]['percentage']
            return f"{sample}\n(H: {percentage:.1f}%)"
        return sample
    if node_id in hybrid_data and hybrid_data[node_id]['percentage']>30:
        percentage = hybrid_data[node_id]['percentage']
        return f"{node_id}\n(H: {percentage:.1f}%)"
    return node_id

def get_node_colors(nodes, hybrids):
    colors=[]
    for node in nodes:
        if node in hybrids:
            colors.append('#ff6b6b')
        else:
            colors.append('#4ecdc4')
    return colors

# === Step 8: Layout and plot ===
plt.figure(figsize=(16,12))
fig, ax = plt.subplots(1,1,figsize=(16,12))
fig.patch.set_facecolor('white')

# Choose layout
if mode == "threshold":
    if layout_choice=="similarity":
        mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
        coords = mds.fit_transform(dm)
        pos = {alignment[i].id:(coords[i,0], coords[i,1]) for i in range(len(alignment)) if alignment[i].id in G.nodes}
    else:
        pos = nx.spring_layout(G, k=0.5, seed=42)
else:  # MSN
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    coords = mds.fit_transform(dm)
    pos = {alignment[i].id:(coords[i,0], coords[i,1]) for i in range(len(alignment)) if alignment[i].id in G.nodes}

# Add hybrid nodes near parents
for hybrid in hybrids:
    if hybrid not in pos:
        parent_positions = [pos[p] for p in G.neighbors(hybrid) if p in pos]
        if parent_positions:
            avg_x = np.mean([p[0] for p in parent_positions])
            avg_y = np.mean([p[1] for p in parent_positions])
            pos[hybrid] = (avg_x + np.random.normal(0,0.05), avg_y + np.random.normal(0,0.05))
        else:
            pos[hybrid] = (np.random.normal(0,1), np.random.normal(0,1))

node_colors = get_node_colors(G.nodes, hybrids)
clean_labels = {node: clean_label(node) for node in G.nodes}

nx.draw_networkx_edges(G,pos,edge_color='lightgray',width=1.5,alpha=0.6)
nx.draw_networkx_nodes(G,pos,node_color=node_colors,node_size=800,edgecolors='black',linewidths=1)
nx.draw_networkx_labels(G,pos,labels=clean_labels,font_size=9,font_weight='bold')

ax.set_title("Cactus Lineage Network Analysis (ITS2+LSU Sequences)",fontsize=16,fontweight='bold')
legend_elements = [Patch(facecolor='#4ecdc4', edgecolor='black', label='Standalone Strains'),
                   Patch(facecolor='#ff6b6b', edgecolor='black', label='Shared Genetics Strains')]
ax.legend(handles=legend_elements, loc='upper right', frameon=True, fancybox=True, shadow=True)
ax.set_axis_off()
plt.tight_layout(pad=3.0)
plt.show()

output_path = os.path.join(BASE_PATH, "Cactus_Lineage_Network_Final.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Network saved as {output_path}")

# === Step 9: Print hybrid summary ===
print("\n=== Hybrid Summary ===")
print(f"Total nodes: {len(G.nodes)}")
print(f"Total edges: {len(G.edges)}")
print(f"Hybrid nodes: {len(hybrids)}")
print(f"Connected components: {nx.number_connected_components(G)}")
