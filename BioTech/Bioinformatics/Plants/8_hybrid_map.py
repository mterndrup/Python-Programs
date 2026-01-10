import os
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# === INPUT FILES ===
BASE_PATH = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus-All"
FASTA_FILE = os.path.join(BASE_PATH, "Cactus_All_top3_ITS_LSU_combined_select.fasta")
#FASTA_FILE = os.path.join(BASE_PATH, "Cactus_All_top3_ITS_LSU_combined_select_lavallo.fasta")
MAFFT_PATH = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Tools\MAFFT\mafft-7.526-win64-signed (1)\mafft-win\mafft.bat"

# Prompt User - Use Most Central Cluster, or All Sequences
choice = input("Show all nodes or one per sample? (all/one): ").strip().lower()
keep_all = choice == "all"

#MSN =
#Threshold =
network_choice = input("Generate Minimum Spanning Network or Threshold Network? (msn/threshold): ").strip().lower()
use_msn = network_choice == "msn"

# Threshold settings
threshold = 0.05
hybrid_threshold = 0.6  # 60% likelihood threshold for calling something a hybrid

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

# === Step 2: Write sequences for MAFFT ===
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


# === Step 5: Distance matrix ===
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

# === Step 6: Build network ===
G = nx.Graph()
for rec in alignment:
    G.add_node(rec.id)

n = len(alignment)

if use_msn:
    print("Building Minimum Spanning Network...")
    for i in range(n):
        for j in range(i + 1, n):
            G.add_edge(alignment[i].id, alignment[j].id, weight=dm[i, j])
    G = nx.minimum_spanning_tree(G, weight='weight')
    print(f"MSN: {len(G.nodes)} nodes, {len(G.edges)} edges")
else:
    print(f"Building Threshold Network (threshold={threshold})...")
    for i in range(n):
        for j in range(i + 1, n):
            if dm[i, j] <= threshold:
                G.add_edge(alignment[i].id, alignment[j].id, weight=1 - dm[i, j])


# === Step 7: Hybrid likelihood ===
def calculate_hybrid_likelihood(seq_idx, distance_matrix, threshold=0.05, connected_nodes=None):
    distances = distance_matrix[seq_idx]
    n_seqs = len(distances)

    # Determine potential parents
    if connected_nodes is not None:  # MSN mode
        potential_parents = [(j, distances[j]) for j in connected_nodes if j != seq_idx]
    else:  # Threshold mode
        potential_parents = [(j, distances[j]) for j in range(n_seqs) if j != seq_idx and 0 < distances[j] <= threshold]

    if len(potential_parents) < 2:
        return 0.0, []

    potential_parents.sort(key=lambda x: x[1])
    likelihood_factors = []

    # Factor 1: number of close relatives
    n_close = len(potential_parents)
    close_factor = min(n_close / 3.0, 1.0)
    likelihood_factors.append(close_factor)

    # Factor 2: distance distribution
    if n_close >= 2:
        distances_to_parents = [d for _, d in potential_parents[:3]]
        avg_distance = np.mean(distances_to_parents)
        distance_variance = np.var(distances_to_parents)
        intermediate_factor = 1.0 - abs(avg_distance - threshold / 2) / (threshold / 2)
        intermediate_factor = max(0, intermediate_factor)
        variance_factor = 1.0 - min(distance_variance / (threshold / 4), 1.0)
        likelihood_factors.extend([intermediate_factor, variance_factor])

    # Factor 3: closest/second-closest ratio
    if len(potential_parents) >= 2:
        ratio = potential_parents[0][1] / potential_parents[1][1] if potential_parents[1][1] > 0 else 0
        likelihood_factors.append(ratio)

    weights = [0.3, 0.25, 0.25, 0.2][:len(likelihood_factors)]
    weights = [w / sum(weights) for w in weights]
    likelihood = min(sum(f * w for f, w in zip(likelihood_factors, weights)), 1.0)

    top_parents = potential_parents[:2]
    return likelihood, top_parents


# Map sequence IDs to indices
id_to_index = {rec.id: i for i, rec in enumerate(alignment)}

# Calculate hybrid likelihoods
hybrid_data = {}
hybrids = []
nodes_to_remove = []

for i in range(n):
    seq_id = alignment[i].id

    # Determine connected nodes if MSN
    connected_nodes = None
    if use_msn:
        neighbors = list(G.neighbors(seq_id))
        connected_nodes = [id_to_index[n] for n in neighbors if n in id_to_index]

    likelihood, parents = calculate_hybrid_likelihood(i, dm, threshold, connected_nodes)

    if likelihood > 0:
        hybrid_data[seq_id] = {'likelihood': likelihood, 'parents': parents, 'percentage': likelihood * 100}

    if likelihood >= hybrid_threshold and len(parents) >= 2:
        sample_name = seq_id.split('_')[0]
        hybrid_name = f"Hybrid_{sample_name}"
        G.add_node(hybrid_name)
        for parent_idx, distance in parents:
            parent_id = alignment[parent_idx].id
            G.add_edge(parent_id, hybrid_name, weight=1 - distance)
        hybrids.append(hybrid_name)
        hybrid_data[hybrid_name] = hybrid_data[seq_id]
        nodes_to_remove.append(seq_id)

for node in nodes_to_remove:
    if G.has_node(node):
        G.remove_node(node)
        print(f"Removed original node {node}, replaced with hybrid version")


# === Helper functions ===
def clean_label(node_id):
    if node_id.startswith("Hybrid_"):
        base_name = node_id.replace("Hybrid_", "")
        if node_id in hybrid_data:
            return f"{base_name}\n({hybrid_data[node_id]['percentage']:.1f}%)"
        return base_name
    parts = node_id.split("_")
    if len(parts) >= 2:
        sample = parts[0]
        if node_id in hybrid_data and hybrid_data[node_id]['percentage'] > 30:
            return f"{sample}\n(H: {hybrid_data[node_id]['percentage']:.1f}%)"
        return sample
    if node_id in hybrid_data and hybrid_data[node_id]['percentage'] > 30:
        return f"{node_id}\n(H: {hybrid_data[node_id]['percentage']:.1f}%)"
    return node_id


def get_node_colors(nodes, hybrids):
    return ['#ff6b6b' if n in hybrids else '#4ecdc4' for n in nodes]


# === Step 8: Visualization ===
plt.style.use('default')
fig, ax = plt.subplots(1, 1, figsize=(16, 12))
fig.patch.set_facecolor('white')

pos = nx.kamada_kawai_layout(G, scale=3) if len(G.nodes) > 20 else nx.spring_layout(G, k=2, iterations=100, seed=42)

clean_labels = {node: clean_label(node) for node in G.nodes}
node_colors = get_node_colors(G.nodes, hybrids)

nx.draw_networkx_edges(G, pos, edge_color='lightgray', width=1.5, alpha=0.6, ax=ax)
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=800, alpha=0.9, edgecolors='black', linewidths=1,
                       ax=ax)

# Smart label positions
y_coords = [pos[n][1] for n in G.nodes]
y_median = np.median(y_coords)
label_pos = {n: (pos[n][0], pos[n][1] - 0.20 if pos[n][1] > y_median else pos[n][1] + 0.15) for n in G.nodes}

nx.draw_networkx_labels(G, label_pos, labels=clean_labels, font_size=9, font_weight='bold', font_family='Arial', ax=ax)

ax.set_title("Cactus Lineage Network Analysis (ITS2+LSU Sequences)", fontsize=16, fontweight='bold', pad=20)
legend_elements = [
    Patch(facecolor='#4ecdc4', edgecolor='black', label='Standalone Strains'),
    Patch(facecolor='#ff6b6b', edgecolor='black', label='Shared Genetics Strains')
]
ax.legend(handles=legend_elements, loc='upper right', frameon=True, fancybox=True, shadow=True)

if not use_msn:
    ax.text(0.02, 0.02, f'Distance Threshold: {threshold}', transform=ax.transAxes,
            fontsize=10, bbox=dict(boxstyle="round,pad=0.3", facecolor='lightblue', alpha=0.7))
ax.text(0.02, 0.08, '% = Shared Genetics Confidence', transform=ax.transAxes,
        fontsize=10, bbox=dict(boxstyle="round,pad=0.3", facecolor='lightyellow', alpha=0.7))
ax.set_axis_off()
plt.tight_layout(pad=3.0)
plt.show()

# Save figure
output_path = os.path.join(BASE_PATH, "Cactus_Lineage_Network_Clean_Fixed.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Network saved as {output_path}")

# === Hybrid analysis table ===
print("\n=== Network Summary ===")
print(
    f"Total nodes: {len(G.nodes)}, edges: {len(G.edges)}, hybrid nodes: {len(hybrids)}, connected components: {nx.number_connected_components(G)}")

if hybrid_data:
    sorted_hybrids = sorted(hybrid_data.items(), key=lambda x: x[1]['percentage'], reverse=True)
    print("\n=== Hybrid Analysis ===")
    print(f"{'Sequence':<20} {'Likelihood':<12} {'Potential Parents':<40} {'Status'}")
    for seq_id, data in sorted_hybrids:
        likelihood_pct = data['percentage']
        status = "HYBRID" if likelihood_pct >= hybrid_threshold * 100 else "Candidate"
        parents_str = ", ".join([f"{alignment[idx].id.split('_')[0]}({dist:.3f})" for idx, dist in data['parents']]) if \
        data['parents'] else "None"
        print(f"{seq_id:<20} {likelihood_pct:>8.1f}%    {parents_str:<40} {status}")
else:
    print("No potential hybrids detected.")
