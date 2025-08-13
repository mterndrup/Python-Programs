import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from adjustText import adjust_text
from sklearn.cluster import DBSCAN
import matplotlib.patches as patches
import re
from scipy.spatial import ConvexHull
import matplotlib.cm as cm

# --- Function definitions ---

def extract_first_element_symbol(col_name):
    match = re.match(r'^([A-Z][a-z]?)', col_name)
    if match:
        return match.group(1)
    return None

k_alpha_energies = {
    'Al': 1.486, 'As': 10.54, 'Au': 9.71, 'Ba': 4.47, 'Bi': 10.84, 'Ca': 3.69,
    'Cd': 23.17, 'Ce': 4.84, 'Cl': 2.622, 'Co': 6.93, 'Cr': 5.41, 'Cu': 8.04,
    'Fe': 6.40, 'Ga': 9.25, 'Hf': 14.90, 'Hg': 9.99, 'K': 3.31, 'La': 4.65,
    'Mg': 1.25, 'Mn': 5.90, 'Mo': 17.48, 'Nb': 16.61, 'Ni': 7.48, 'P': 2.01,
    'Pb': 10.55, 'Pd': 21.18, 'Pt': 9.44, 'Rb': 13.39, 'S': 2.31, 'Sb': 26.39,
    'Se': 11.22, 'Si': 1.74, 'Sn': 25.27, 'Sr': 14.16, 'Ta': 57.54, 'Te': 27.47,
    'Th': 59.31, 'Ti': 4.51, 'Tl': 10.48, 'U': 98.43, 'V': 4.95, 'W': 59.32,
    'Y': 14.96, 'Zn': 8.64, 'Zr': 15.77,
}

def get_k_alpha_energy_keV(symbol):
    return k_alpha_energies.get(symbol)

# --- Load and preprocess data ---
file_path = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Samples_NoErr.csv"
df = pd.read_csv(file_path)

df['File #'] = df['File #'].astype(str).str.strip()

# Remove samples #357 and #366 only (keep #13)
df = df[~df['File #'].isin(['357', '366'])].reset_index(drop=True)

df.replace(["< LOD", ""], np.nan, inplace=True)

exclude_cols = ['File #', 'Location']
concentration_cols = [col for col in df.columns if col not in exclude_cols]

# Convert all concentration columns to numeric
for col in concentration_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

df_filled = df[concentration_cols].fillna(0)

# No filtering step here — keep all elements as-is

# --- Z-score normalization ---
scaler = StandardScaler()
scaled_data = scaler.fit_transform(df_filled)

# --- PCA ---
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_data)
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['SampleNum'] = range(1, len(pca_df) + 1)
pca_df['Location'] = df['Location']

# --- Plot 1: PCA Scatter Plot ---
plt.figure(figsize=(12, 6))

locations = pca_df['Location'].unique()
colors_map = {loc: cm.tab10(i % 10) for i, loc in enumerate(locations)}
colors = pca_df['Location'].map(colors_map)

plt.scatter(pca_df['PC1'], pca_df['PC2'], c=colors, alpha=0.8, s=60, edgecolor='k', linewidth=0.5, zorder=2)

texts = []
for i, num in enumerate(pca_df['SampleNum']):
    texts.append(plt.text(pca_df['PC1'][i], pca_df['PC2'][i], str(num), fontsize=8, fontweight='bold', zorder=3))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', lw=0.5))

plt.title('PCA of XRF Samples')
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% variance)")
plt.grid(True, zorder=1)

handles = [plt.Line2D([0], [0], marker='o', color='w', label=f"{i+1}: {loc}",
                      markerfacecolor=colors_map[loc], markersize=8) for i, loc in enumerate(locations)]
plt.legend(handles=handles, title='Location (Label Number: Location)', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

# --- Plot 2: PCA Loadings with Cluster Outlines ---

loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=concentration_cols)

fig, ax = plt.subplots(figsize=(12, 8))
colors = []
for element in loadings.index:
    val = loadings.loc[element, 'PC1']
    if val > 0.5:
        colors.append('red')
    elif val < -0.5:
        colors.append('blue')
    else:
        colors.append('black')

ax.scatter(loadings['PC1'], loadings['PC2'], color=colors, zorder=2)
texts = []
for element in loadings.index:
    texts.append(ax.text(loadings.loc[element, 'PC1'], loadings.loc[element, 'PC2'], element,
                         fontsize=9, fontweight='bold', color='black', zorder=3))
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=0.5, color='gray'))

X = loadings[['PC1', 'PC2']].values
dbscan = DBSCAN(eps=0.3, min_samples=2)
clusters = dbscan.fit_predict(X)
unique_clusters = set(clusters)
cmap = plt.cm.get_cmap('tab10', max(len(unique_clusters), 10))

for cluster_id in unique_clusters:
    if cluster_id == -1:
        continue
    cluster_points = X[clusters == cluster_id]
    if len(cluster_points) < 3:
        continue
    hull = ConvexHull(cluster_points)
    hull_points = cluster_points[hull.vertices]
    polygon = patches.Polygon(hull_points, closed=True, facecolor=cmap(cluster_id % 10),
                              alpha=0.2, edgecolor=cmap(cluster_id % 10), linewidth=2, zorder=1)
    ax.add_patch(polygon)

ax.axhline(0, color='black', linewidth=1, zorder=1)
ax.axvline(0, color='black', linewidth=1, zorder=1)
ax.set_xlabel('Major Loading Trend (PC1)')
ax.set_ylabel('Orthogonal Loading Trend (PC2)')
ax.set_title('Element Correlations from LA County Soils')
ax.grid(True, zorder=0)
plt.tight_layout()
plt.show()

# --- Plot 3: Combined Metric (Presence Ratio × Mean Concentration) vs Energy (keV) ---

presence_ratio = (df_filled > 0).sum() / len(df_filled)
mean_ppm = df_filled.mean()
combined_metric = presence_ratio * mean_ppm

element_keV = {}
filtered_elements = []

for col in combined_metric.index:
    symbol = extract_first_element_symbol(col)
    energy = get_k_alpha_energy_keV(symbol) if symbol else None
    if energy is not None:
        element_keV[col] = energy
        filtered_elements.append(col)
    else:
        print(f"Warning: No K-alpha energy for element '{symbol}' (column '{col}') - skipping")

sorted_elements = sorted(element_keV.items(), key=lambda x: x[1])

energies = [e for _, e in sorted_elements]
metrics = [combined_metric[el] for el, _ in sorted_elements]

plt.figure(figsize=(14, 6))
plt.scatter(energies, metrics, color='purple', s=100, alpha=0.8)

for el, x, y in zip([el for el, _ in sorted_elements], energies, metrics):
    plt.text(x, y * 1.05, el, ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.xlabel('Energy (keV)')
plt.ylabel('Presence Ratio × Mean Concentration (ppm)')
plt.title('Element Presence and Abundance vs Emission Energy')
plt.grid(True)

energy_min = min(energies) * 0.8
energy_max = max(energies) * 1.2
plt.xlim(energy_min, energy_max)
plt.tight_layout()
plt.show()

# --- Plot 4: Correlate matrix elements with combined metric from Plot 3 ---

matrix_elements = ['SiO2', 'Fe', 'Ca', 'K2O', 'Al2O3']
matrix_elements = [elem for elem in matrix_elements if elem in df_filled.columns]

matrix_means = df_filled[matrix_elements].mean()

oxide_to_element = {
    'SiO2': 'Si',
    'Fe': 'Fe',
    'Ca': 'Ca',
    'K2O': 'K',
    'Al2O3': 'Al'
}

matrix_combined_metric = {}
matrix_energies = {}

for oxide in matrix_elements:
    element = oxide_to_element.get(oxide)
    if element and element in combined_metric.index:
        matrix_combined_metric[oxide] = combined_metric[element]
        energy = get_k_alpha_energy_keV(element)
        if energy:
            matrix_energies[oxide] = energy

fig, ax = plt.subplots(figsize=(10, 6))
x = [matrix_means[oxide] for oxide in matrix_elements if oxide in matrix_means.index]
y = [matrix_combined_metric.get(oxide, np.nan) for oxide in matrix_elements if oxide in matrix_means.index]
colors = [matrix_energies.get(oxide, 0) for oxide in matrix_elements if oxide in matrix_means.index]

sc = ax.scatter(x, y, c=colors, cmap='viridis', s=150, edgecolor='k', alpha=0.8)
for i, oxide in enumerate(matrix_elements):
    if oxide in matrix_means.index:
        ax.text(x[i], y[i]*1.05, oxide, ha='center', fontsize=12, fontweight='bold')

cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('K-alpha Energy (keV)')

ax.set_xlabel('Mean Concentration of Matrix Element (ppm)')
ax.set_ylabel('Presence Ratio × Mean Concentration (ppm)')
ax.set_title('Correlation of Matrix Element Concentration to Combined Metric with Energy')

ax.grid(True)
plt.tight_layout()
plt.show()
