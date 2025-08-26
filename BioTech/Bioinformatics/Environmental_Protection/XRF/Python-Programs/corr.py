import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Path to your CSV file
file_path = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Samples_NoErr.csv"

# Load CSV, replace '< LOD' with NaN
df = pd.read_csv(file_path)
df.replace("< LOD", np.nan, inplace=True)

# Select element concentration columns (exclude 'File #' and 'Location')
concentration_cols = [col for col in df.columns if col not in ['File #', 'Location']]

# Convert all concentration columns to numeric
for col in concentration_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Calculate Pearson and Spearman correlation matrices
pearson_corr = df[concentration_cols].corr(method='pearson')
spearman_corr = df[concentration_cols].corr(method='spearman')

print("Pearson Correlation Matrix (Concentrations):")
print(pearson_corr.round(2))

print("\nSpearman Correlation Matrix (Concentrations):")
print(spearman_corr.round(2))

# Plot heatmaps with all elements clearly visible
plt.figure(figsize=(20, 18))  # Larger figure to fit all labels

# Pearson
plt.subplot(1, 2, 1)
sns.heatmap(pearson_corr, annot=False, cmap='coolwarm', center=0,
            xticklabels=pearson_corr.columns, yticklabels=pearson_corr.columns)
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.title('Pearson Correlation', fontsize=14)

# Spearman
plt.subplot(1, 2, 2)
sns.heatmap(spearman_corr, annot=False, cmap='coolwarm', center=0,
            xticklabels=spearman_corr.columns, yticklabels=spearman_corr.columns)
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.title('Spearman Correlation', fontsize=14)

plt.tight_layout()
plt.show()
