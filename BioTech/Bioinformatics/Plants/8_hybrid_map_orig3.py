import os
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import seaborn as sns

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
        for j in range(i + 1, n):
            seq1 = str(aligned_records[i].seq)
            seq2 = str(aligned_records[j].seq)
            matches = sum(a == b for a, b in zip(seq1, seq2))
            distance = 1 - matches / len(seq1)
            dm[i, j] = dm[j, i] = distance
    return dm


dm = build_distance_matrix(alignment)

# === Step 6: Build network graph ===
G = nx.Graph()  # Use undirected graph for cleaner layout
for rec in alignment:
    G.add_node(rec.id)

# Connect sequences below threshold
threshold = 0.05  # adjust as needed
hybrid_threshold = 0.6  # 60% likelihood threshold for calling something a hybrid
n = len(alignment)
for i in range(n):
    for j in range(i + 1, n):
        if dm[i, j] <= threshold:
            G.add_edge(alignment[i].id, alignment[j].id, weight=1 - dm[i, j])


# === Step 7: Calculate hybrid likelihood and add hybrid nodes ===
def calculate_hybrid_likelihood(seq_idx, distance_matrix, threshold=0.05):
    """Calculate likelihood that a sequence is a hybrid based on distance patterns"""
    distances = distance_matrix[seq_idx]
    n_seqs = len(distances)

    # Find potential parents (sequences within threshold, excluding self)
    potential_parents = [(j, distances[j]) for j in range(n_seqs)
                         if j != seq_idx and 0 < distances[j] <= threshold]

    if len(potential_parents) < 2:
        return 0.0, []

    # Sort by distance (closest first)
    potential_parents.sort(key=lambda x: x[1])

    # Calculate hybrid likelihood based on multiple factors
    likelihood_factors = []

    # Factor 1: Number of close relatives (more = higher hybrid likelihood)
    n_close = len(potential_parents)
    close_factor = min(n_close / 3.0, 1.0)  # Normalize to max of 1.0
    likelihood_factors.append(close_factor)

    # Factor 2: Distance distribution pattern
    # Hybrids often show intermediate distances to multiple parents
    if n_close >= 2:
        distances_to_parents = [d for _, d in potential_parents[:3]]  # Top 3 closest
        avg_distance = np.mean(distances_to_parents)
        distance_variance = np.var(distances_to_parents)

        # Intermediate distances (not too close, not too far) suggest hybridization
        intermediate_factor = 1.0 - abs(avg_distance - threshold / 2) / (threshold / 2)
        intermediate_factor = max(0, intermediate_factor)

        # Low variance in distances to multiple parents suggests hybrid
        variance_factor = 1.0 - min(distance_variance / (threshold / 4), 1.0)

        likelihood_factors.extend([intermediate_factor, variance_factor])

    # Factor 3: Genetic distance ratio test
    # Compare distance to closest vs second closest (hybrids show similar distances)
    if len(potential_parents) >= 2:
        closest_dist = potential_parents[0][1]
        second_closest_dist = potential_parents[1][1]
        ratio = closest_dist / second_closest_dist if second_closest_dist > 0 else 0
        ratio_factor = ratio  # Higher ratio = more similar distances = more likely hybrid
        likelihood_factors.append(ratio_factor)

    # Calculate overall likelihood (weighted average)
    weights = [0.3, 0.25, 0.25, 0.2]  # Adjust weights as needed
    weights = weights[:len(likelihood_factors)]
    weights = [w / sum(weights) for w in weights]  # Normalize weights

    likelihood = sum(f * w for f, w in zip(likelihood_factors, weights))
    likelihood = min(likelihood, 1.0)  # Cap at 100%

    # Return top 2 potential parents
    top_parents = potential_parents[:2]

    return likelihood, top_parents


# Calculate hybrid likelihoods for all sequences
hybrid_data = {}
hybrids = []
nodes_to_remove = []  # Track original nodes that become hybrids

for i in range(n):
    likelihood, parents = calculate_hybrid_likelihood(i, dm, threshold)
    seq_id = alignment[i].id

    if likelihood > 0:
        hybrid_data[seq_id] = {
            'likelihood': likelihood,
            'parents': parents,
            'percentage': likelihood * 100
        }

    # Add to network if likelihood is high enough
    if likelihood >= hybrid_threshold and len(parents) >= 2:
        sample_name = seq_id.split('_')[0]
        hybrid_name = f"Hybrid_{sample_name}"
        G.add_node(hybrid_name)

        # Add edges to potential parents
        for parent_idx, distance in parents:
            parent_id = alignment[parent_idx].id
            G.add_edge(parent_id, hybrid_name, weight=1 - distance)

        hybrids.append(hybrid_name)

        # Also store hybrid data for the hybrid node
        hybrid_data[hybrid_name] = hybrid_data[seq_id]

        # Mark original node for removal
        nodes_to_remove.append(seq_id)

# Remove original nodes that now have hybrid representations
for node in nodes_to_remove:
    if G.has_node(node):
        G.remove_node(node)
        print(f"Removed original node {node}, replaced with hybrid version")


# === Helper functions for cleaner labels (FIXED VERSION) ===
def clean_label(node_id):
    """Clean up node labels for better readability - Fixed to handle all cases"""
    if node_id.startswith("Hybrid_"):
        base_name = node_id.replace("Hybrid_", "")
        # Add likelihood percentage if available
        if node_id in hybrid_data:
            percentage = hybrid_data[node_id]['percentage']
            return f"{base_name}\n({percentage:.1f}%)"
        return base_name

    # Extract key parts from sequence ID
    parts = node_id.split("_")
    if len(parts) >= 2:
        sample = parts[0]
        # Add hybrid likelihood for regular sequences if they have significant likelihood
        if node_id in hybrid_data and hybrid_data[node_id]['percentage'] > 30:
            percentage = hybrid_data[node_id]['percentage']
            return f"{sample}\n(H: {percentage:.1f}%)"
        return sample

    # Fallback for any unexpected node ID format
    if node_id in hybrid_data and hybrid_data[node_id]['percentage'] > 30:
        percentage = hybrid_data[node_id]['percentage']
        return f"{node_id}\n(H: {percentage:.1f}%)"

    return node_id  # Return the full ID if nothing else works


def get_node_colors(nodes, hybrids):
    """Generate distinct colors for different node types"""
    colors = []
    for node in nodes:
        if node in hybrids:
            colors.append('#ff6b6b')  # Red for hybrids
        else:
            # Use different shades of blue for regular nodes
            colors.append('#4ecdc4')  # Teal for regular nodes
    return colors


# === Step 8: Create improved visualization (FIXED VERSION) ===
plt.style.use('default')
fig, ax = plt.subplots(1, 1, figsize=(16, 12))
fig.patch.set_facecolor('white')

# Use better layout algorithms
if len(G.nodes) > 20:
    pos = nx.kamada_kawai_layout(G, scale=3)
else:
    pos = nx.spring_layout(G, k=2, iterations=100, seed=42)

# Clean up node labels - ENSURE ALL NODES GET LABELS
clean_labels = {}
for node in G.nodes:
    clean_labels[node] = clean_label(node)
    # Debug: print to verify all nodes have labels
    print(f"Node: {node} -> Label: {clean_labels[node]}")

# Verify we have labels for all nodes
print(f"Total nodes: {len(G.nodes)}")
print(f"Total labels: {len(clean_labels)}")
print(f"Node list: {list(G.nodes)}")

# Get node colors
node_colors = get_node_colors(G.nodes, hybrids)

# Draw edges first (behind nodes)
nx.draw_networkx_edges(G, pos,
                       edge_color='lightgray',
                       width=1.5,
                       alpha=0.6,
                       ax=ax)

# Draw nodes with better styling
nx.draw_networkx_nodes(G, pos,
                       node_color=node_colors,
                       node_size=800,
                       alpha=0.9,
                       edgecolors='black',
                       linewidths=1,
                       ax=ax)

# Draw labels with smart positioning - above or below based on node position
label_pos = {}
# Get the y-coordinate range to determine what's "top" vs "bottom"
y_coords = [pos[node][1] for node in pos.keys()]
y_median = np.median(y_coords)

for node in G.nodes:
    if node in pos:
        x, y = pos[node]
        # If node is in top half of plot, put label below; otherwise above
        if y > y_median:
            label_pos[node] = (x, y - 0.20)  # Below for top nodes
        else:
            label_pos[node] = (x, y + 0.15)  # Above for bottom nodes
    else:
        print(f"Warning: Node {node} not found in position dictionary")

print(f"Drawing labels for {len(clean_labels)} nodes with smart positioning")

# Draw all labels with smart positioning
nx.draw_networkx_labels(G, label_pos,
                        labels=clean_labels,
                        font_size=9,
                        font_weight='bold',
                        font_family='Arial',
                        ax=ax)

# Add title and styling
ax.set_title("Cactus Lineage Network Analysis (ITS2+LSU Sequences)",
             fontsize=16,
             fontweight='bold',
             pad=20)

# Create legend
from matplotlib.patches import Patch

legend_elements = [
    Patch(facecolor='#4ecdc4', edgecolor='black', label='Standalone Strains'),
    Patch(facecolor='#ff6b6b', edgecolor='black', label='Shared Genetics Strains')
]
ax.legend(handles=legend_elements, loc='upper right', frameon=True, fancybox=True, shadow=True)

# Add distance threshold info
ax.text(0.02, 0.02, f'Distance Threshold: {threshold}',
        transform=ax.transAxes,
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor='lightblue', alpha=0.7))

# Add percentage explanation
ax.text(0.02, 0.08, '% = Shared Genetics Confidence',
        transform=ax.transAxes,
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor='lightyellow', alpha=0.7))

# Remove axes
ax.set_axis_off()

# Adjust layout with more padding for title
plt.tight_layout(pad=3.0)

# Display the plot
plt.show()

# Save with higher DPI
output_path = os.path.join(BASE_PATH, "Cactus_Lineage_Network_Clean_Fixed.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Fixed network saved as {output_path}")

# === Optional: Create detailed hybrid analysis table ===
print("\n=== Network Summary ===")
print(f"Total nodes: {len(G.nodes)}")
print(f"Total edges: {len(G.edges)}")
print(f"Hybrid nodes: {len(hybrids)}")
print(f"Connected components: {nx.number_connected_components(G)}")

# Create comprehensive hybrid analysis
print(f"\n=== Hybrid Analysis (Threshold: {hybrid_threshold * 100}% likelihood) ===")
if hybrid_data:
    # Sort by likelihood
    sorted_hybrids = sorted(hybrid_data.items(), key=lambda x: x[1]['percentage'], reverse=True)

    print(f"{'Sequence':<20} {'Likelihood':<12} {'Potential Parents':<40} {'Status'}")
    print("-" * 90)

    for seq_id, data in sorted_hybrids:
        likelihood_pct = data['percentage']
        status = "HYBRID" if likelihood_pct >= hybrid_threshold * 100 else "Candidate"

        # Format parent information
        if data['parents']:
            parent_info = []
            for parent_idx, distance in data['parents']:
                if seq_id.startswith("Hybrid_"):
                    parent_name = alignment[parent_idx].id.split('_')[0]
                else:
                    parent_name = alignment[parent_idx].id.split('_')[0]
                parent_info.append(f"{parent_name}({distance:.3f})")
            parents_str = ", ".join(parent_info)
        else:
            parents_str = "None identified"

        seq_display = seq_id.split('_')[0] if not seq_id.startswith("Hybrid_") else seq_id
        print(f"{seq_display:<20} {likelihood_pct:>8.1f}%    {parents_str:<40} {status}")

    # Summary statistics
    high_likelihood = [d for d in hybrid_data.values() if d['percentage'] >= 70]
    medium_likelihood = [d for d in hybrid_data.values() if 40 <= d['percentage'] < 70]
    low_likelihood = [d for d in hybrid_data.values() if d['percentage'] < 40]

    print(f"\n=== Hybrid Likelihood Distribution ===")
    print(f"High likelihood (≥70%): {len(high_likelihood)} sequences")
    print(f"Medium likelihood (40-69%): {len(medium_likelihood)} sequences")
    print(f"Low likelihood (<40%): {len(low_likelihood)} sequences")
    print(f"Average likelihood: {np.mean([d['percentage'] for d in hybrid_data.values()]):.1f}%")

else:
    print("No potential hybrids detected with current parameters.")

print(f"\n=== Parameters Used ===")
print(f"Distance threshold: {threshold}")
print(f"Hybrid likelihood threshold: {hybrid_threshold * 100}%")
print(f"Minimum parents required: 2")