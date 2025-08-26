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

# --- Ask user which chart to generate ---
print("Select a chart to generate:")
print("1: PCA Scatter Plot of Samples")
print("2: PCA Loadings with Clusters")
print("3: Presence Ratio × Mean Concentration vs Energy")
print("4: Correlation of Matrix Elements with Combined Metric")
print("0: Generate All Charts")

chart_choice = input("Enter the number of the chart to generate (0-4): ").strip()

# --- Load and preprocess data ---
file_path = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Samples_NoErr.csv"
# --- Load full dataset ---
df_all = pd.read_csv(file_path)
df_all['File #'] = df_all['File #'].astype(str).str.strip()
df_all.replace(["< LOD", ""], np.nan, inplace=True)

# Concentration columns
exclude_cols = ['File #', 'Location']
concentration_cols_all = [col for col in df_all.columns if col not in exclude_cols]

# Convert to numeric
for col in concentration_cols_all:
    df_all[col] = pd.to_numeric(df_all[col], errors='coerce')

# Fill missing values and log-transform
df_filled_all = df_all[concentration_cols_all].fillna(0)
df_log_all = np.log10(df_filled_all + 1e-6)

# Z-score normalization for Chart 2
scaler_all = StandardScaler()
scaled_data_all = scaler_all.fit_transform(df_log_all)

# PCA for Chart 2
pca_all = PCA(n_components=2)
pca_result_all = pca_all.fit_transform(scaled_data_all)
pca_df_all = pd.DataFrame(pca_result_all, columns=['PC1', 'PC2'])
pca_df_all['SampleNum'] = range(1, len(pca_df_all) + 1)
pca_df_all['Location'] = df_all['Location']

# --- Chart 1 (filtered dataset) ---
df = df_all[~df_all['File #'].isin(['357', '366'])].reset_index(drop=True)
concentration_cols = [col for col in df.columns if col not in exclude_cols]
df_filled = df[concentration_cols].fillna(0)
df_log = np.log10(df_filled + 1e-6)

scaler = StandardScaler()
scaled_data = scaler.fit_transform(df_log)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_data)
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['SampleNum'] = range(1, len(pca_df) + 1)
pca_df['Location'] = df['Location']


if chart_choice not in ['0', '1', '2', '3', '4']:
    print("⚠️ Invalid choice, defaulting to All Charts")
    chart_choice = '0'

# Example usage later in your plotting code:
if chart_choice in ['0', '1']:
    # --- Plot 1: PCA Scatter Plot ---
    plt.figure(figsize=(12, 6))

    locations = pca_df['Location'].unique()
    colors_map = {loc: cm.tab10(i % 10) for i, loc in enumerate(locations)}
    colors = pca_df['Location'].map(colors_map)

    plt.scatter(pca_df['PC1'], pca_df['PC2'], c=colors, alpha=0.8, s=380, edgecolor='k', linewidth=0.5, zorder=2)

    texts = []
    for i, num in enumerate(pca_df['SampleNum']):
        x = pca_df['PC1'][i]
        y = pca_df['PC2'][i]
        texts.append(plt.text(x, y, str(num),
                              fontsize=10,
                              fontweight='bold',
                              ha='center',  # horizontal alignment
                              va='center',  # vertical alignment
                              color='black',  # optional, for contrast inside marker
                              zorder=3))

    plt.title('PCA of XRF Samples')
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)")
    plt.grid(True, zorder=1)

    handles = [plt.Line2D([0], [0], marker='o', color='w', label=f"{i + 1}: {loc}",
                          markerfacecolor=colors_map[loc], markersize=8) for i, loc in enumerate(locations)]
    plt.legend(handles=handles, title='Location (Label Number: Location)', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()

if chart_choice in ['0', '2']:
    # --- Plot 2: PCA Loadings with Metal vs Natural Colors ---
    loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=concentration_cols)

    # Filter out elements near 0,0
    threshold = 0.05
    filtered_loadings = loadings[(abs(loadings['PC1']) > threshold) | (abs(loadings['PC2']) > threshold)]

    fig, ax = plt.subplots(figsize=(17, 8))

    # Define metal vs natural
    metal_elements = {'Pd', 'Pt', 'Cu', 'Fe', 'Cr', 'Ti', 'Ni', 'Pb', 'Mn', 'Zn', 'Co', 'V', 'Mo', 'Sb'}
    natural_elements = set(filtered_loadings.index) - metal_elements

    # Assign colors
    element_colors = {}
    for el in filtered_loadings.index:
        if el in metal_elements:
            element_colors[el] = 'red'
        else:
            element_colors[el] = 'green'

    # Scatter plot
    for element in filtered_loadings.index:
        ax.scatter(filtered_loadings.loc[element, 'PC1'],
                   filtered_loadings.loc[element, 'PC2'],
                   color=element_colors[element], s=380, zorder=2)

    # Add text labels in the middle of each point
    for element in filtered_loadings.index:
        ax.text(filtered_loadings.loc[element, 'PC1'],
                filtered_loadings.loc[element, 'PC2'],
                element,
                fontsize=9, fontweight='bold', ha='center', va='center',
                color='black', zorder=3)

    # Axes lines, labels, title, grid
    ax.axhline(0, color='black', linewidth=1, zorder=1)
    ax.axvline(0, color='black', linewidth=1, zorder=1)
    ax.set_xlabel('Major Loading Trend (PC1)')
    ax.set_ylabel('Orthogonal Loading Trend (PC2)')
    ax.set_title('Element Correlations from LA County Soils')
    ax.grid(True, zorder=0)

    # --- Legend with filtered elements on the right ---
    from matplotlib.lines import Line2D

    legend_elements = [Line2D([0], [0], marker='o', color='w', label=el,
                              markerfacecolor=element_colors[el], markersize=12)
                       for el in filtered_loadings.index]

    ax.legend(handles=legend_elements, title='Elements', bbox_to_anchor=(1.05, 1),
              loc='upper left', fontsize=9)

    plt.tight_layout()
    plt.show()

if chart_choice in ['0', '3']:
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

if chart_choice in ['0', '4']:
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
            ax.text(x[i], y[i] * 1.05, oxide, ha='center', fontsize=12, fontweight='bold')

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('K-alpha Energy (keV)')

    ax.set_xlabel('Mean Concentration of Matrix Element (ppm)')
    ax.set_ylabel('Presence Ratio × Mean Concentration (ppm)')
    ax.set_title('Correlation of Matrix Element Concentration to Combined Metric with Energy')

    ax.grid(True)
    plt.tight_layout()
    plt.show()



