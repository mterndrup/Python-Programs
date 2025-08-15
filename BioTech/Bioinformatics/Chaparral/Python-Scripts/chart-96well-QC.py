import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# -----------------------------
# User input
# -----------------------------
k = input("Enter k-mer value (6 or 8): ").strip()
if k not in ["6", "8"]:
    raise ValueError("k-mer value must be 6 or 8")

# -----------------------------
# Paths
# -----------------------------
combined_file = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Checkpoints\combined_k{}_extracted.csv".format(k)
charts_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Charts"
checkpoint_file = os.path.join(charts_dir, f"barcode_cluster_map_k{k}.csv")
os.makedirs(charts_dir, exist_ok=True)

# -----------------------------
# Plate layout (H at top, A at bottom)
# -----------------------------
plate_rows = ["H","G","F","E","D","C","B","A"]
plate_cols = list(range(1,13))
plate_data = np.full((8,12), "", dtype=object)

# -----------------------------
# Build plate mapping barcode -> position
# -----------------------------
barcode_counter = 1
barcode_positions = {}  # barcode_num -> (row_idx, col_idx)
for col_idx, col_num in enumerate(plate_cols):
    for row_idx in reversed(range(8)):  # start from bottom A
        barcode_num = f"{barcode_counter:02d}"
        barcode_positions[barcode_num] = (row_idx, col_idx)
        plate_data[row_idx, col_idx] = barcode_num
        barcode_counter += 1

# -----------------------------
# Read combined clusters
# -----------------------------
combined_df = pd.read_csv(combined_file)

# -----------------------------
# Extract barcode from source_file
# -----------------------------
combined_df['barcode'] = combined_df['source_file'].str.extract(r'(barcode\d{2})')[0]

# -----------------------------
# Compute cluster counts per barcode
# -----------------------------
if os.path.exists(checkpoint_file):
    print("Loading checkpoint from:", checkpoint_file)
    df_chk = pd.read_csv(checkpoint_file, index_col=0)
    barcode_fraction_map = {b: {'native_count': row['native_count'], 'contaminant': row['contaminant']}
                            for b,row in df_chk.iterrows() if b in barcode_positions}
else:
    # -----------------------------
    # Determine majority barcode per cluster
    # -----------------------------
    cluster_counts = combined_df.groupby(['cluster','barcode']).size().reset_index(name='count')
    majority_barcode = cluster_counts.loc[cluster_counts.groupby('cluster')['count'].idxmax()]
    majority_barcode_map = dict(zip(majority_barcode['cluster'], majority_barcode['barcode']))

    # -----------------------------
    # Assign native/contaminant per cluster
    # -----------------------------
    cluster_summary = combined_df[['cluster','barcode']].drop_duplicates()
    cluster_summary['is_native'] = cluster_summary.apply(
        lambda row: row['barcode'] == majority_barcode_map.get(row['cluster'], None), axis=1
    )

    # Aggregate **cluster counts** per barcode (key fix)
    barcode_summary = cluster_summary.groupby('barcode')['is_native'].agg(['sum','count'])
    barcode_summary['contaminant'] = barcode_summary['count'] - barcode_summary['sum']
    barcode_summary = barcode_summary.rename(columns={'sum':'native_count'})

    # Only keep barcodes that exist in plate layout
    barcode_fraction_map = {b: {'native_count': row['native_count'], 'contaminant': row['contaminant']}
                            for b,row in barcode_summary.iterrows() if b in barcode_positions}

    # Save checkpoint
    barcode_summary.to_csv(checkpoint_file)
    print("Checkpoint saved to:", checkpoint_file)

# -----------------------------
# Debug: print barcodes with data
# -----------------------------
print("Barcodes with data:", list(barcode_fraction_map.keys()))

# -----------------------------
# Plot 96-well plate with native/contaminant cluster segments
# -----------------------------
fig, ax = plt.subplots(figsize=(12,8))
ax.set_xlim(-0.5,11.5)
ax.set_ylim(-0.5,7.5)

for barcode_num, (row_idx, col_idx) in barcode_positions.items():
    x0, y0 = col_idx-0.5, row_idx-0.5
    x1, y1 = col_idx+0.5, row_idx+0.5
    # Draw outer square
    ax.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0], color="black")
    # Barcode label
    ax.text(col_idx, y1-0.1, barcode_num, ha="center", va="top", fontsize=10)

    # Draw stacked cluster segments
    if barcode_num in barcode_fraction_map:
        native_count = int(barcode_fraction_map[barcode_num]['native_count'])
        contaminant_count = int(barcode_fraction_map[barcode_num]['contaminant'])
        total = native_count + contaminant_count
        if total > 0:
            segment_height = (y1 - y0) / total
            bottom = y0
            # Contaminant (red) at bottom
            for _ in range(contaminant_count):
                ax.fill_between([x0,x1],[bottom,bottom],[bottom+segment_height,bottom+segment_height],
                                color="red", edgecolor="black")
                bottom += segment_height
            # Native (green) on top
            for _ in range(native_count):
                ax.fill_between([x0,x1],[bottom,bottom],[bottom+segment_height,bottom+segment_height],
                                color="green", edgecolor="black")
                bottom += segment_height

# -----------------------------
# Format axes
# -----------------------------
ax.set_xticks(np.arange(12))
ax.set_yticks(np.arange(8))
ax.set_xticklabels(range(1,13))
ax.set_yticklabels(plate_rows)
ax.tick_params(axis='both', top=True, labeltop=True, bottom=False, labelbottom=False)
ax.set_xlabel("Column")
ax.set_ylabel("Row")
ax.set_title(f"96-Well Plate Native vs Contaminant Clusters (k={k})", pad=20)

# -----------------------------
# Save chart
# -----------------------------
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_file = os.path.join(charts_dir, f"96well_plate_k{k}_clusters_{timestamp}.png")
plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.show()
print(f"Chart saved to: {output_file}")
