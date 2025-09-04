from matplotlib_venn import venn3
import matplotlib.pyplot as plt

# Counts from your dataset
counts_summary = {
    "Urban only": 14,
    "Recreational only": 8,
    "Fire only": 0,
    "Urban ∩ Fire": 11,
    "Rec ∩ Fire": 21,
    "Urban ∩ Rec": 0,
    "Urban ∩ Rec ∩ Fire": 0,
}

# Create Venn diagram
plt.figure(figsize=(8, 8))
venn = venn3(
    subsets=(
        counts_summary["Urban only"],
        counts_summary["Recreational only"],
        counts_summary["Urban ∩ Rec"],
        counts_summary["Fire only"],
        counts_summary["Urban ∩ Fire"],
        counts_summary["Rec ∩ Fire"],
        counts_summary["Urban ∩ Rec ∩ Fire"]
    ),
    set_labels=("Urban", "Fire-affected", "Recreational")
)

# Custom colors: Urban = purple, Recreational = green, Fire = red
colors = {
    '100': '#A020F0',  # Urban only (purple)
    '010': '#00CC66',  # Recreational only (green)
    '001': '#FF3333',  # Fire only (red)
    '110': '#9999FF',
    '101': '#FF99CC',
    '011': '#99FF99',
    '111': '#CCCCCC'
}

for subset, color in colors.items():
    patch = venn.get_patch_by_id(subset)
    if patch:
        patch.set_color(color)
        patch.set_alpha(0.6)
        patch.set_edgecolor("black")   # outlines
        patch.set_linewidth(2)

# Remove labels for zero subsets
for subset in ['100', '010', '110', '001', '101', '011', '111']:
    label = venn.get_label_by_id(subset)
    if label and label.get_text() == '0':
        label.set_text('')
    elif label:  # increase font size & make bold
        label.set_fontsize(14)
        label.set_fontweight("bold")

# Adjust bottom labels (set_labels) to align evenly
venn_labels = venn.set_labels
positions = [(-0.6, -0.5), (0, -0.5), (0.6, -0.5)]

for label, pos in zip(venn_labels, positions):
    if label:
        label.set_position(pos)
        label.set_verticalalignment('center')
        label.set_horizontalalignment('center')
        label.set_fontsize(14)
        label.set_fontweight("bold")

# Update title
plt.title("Venn Diagram of All Soil Samples", fontsize=16, fontweight="bold")

# Save as SVG
plt.savefig("soil_samples_venn.svg", format="svg")

plt.show()
