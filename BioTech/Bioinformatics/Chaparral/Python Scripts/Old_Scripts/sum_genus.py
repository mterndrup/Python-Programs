import pandas as pd
import os

# Define base directory
base_dir = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral"

# Define file paths
input_file = os.path.join(base_dir, "dna_counts_final_sorted.csv")
output_file = os.path.join(base_dir, "filtered_genus_b_counts.csv")

# Read the comma-delimited data (default)
df = pd.read_csv(input_file)

# Normalize column names: lowercase, strip spaces
df.columns = df.columns.str.strip().str.lower()

print("Detected columns:", df.columns.tolist())

# Drop rows with missing genus or species
df = df.dropna(subset=["genus", "species"])

# Filter out rows with 'uncultured' or 'sp.' in genus or species
def is_valid(row):
    genus = str(row["genus"]).strip().lower()
    species = str(row["species"]).strip().lower()
    return not (
        "uncultured" in genus or "uncultured" in species or
        species == "sp." or genus == "sp."
    )

filtered_df = df[df.apply(is_valid, axis=1)]

# Define the b* columns to keep
b_columns = ["b05", "b08", "b21", "b22", "b23", "b24", "b62", "b63"]

# Keep only b_columns that actually exist in the DataFrame
existing_b_columns = [col for col in b_columns if col in filtered_df.columns]

# Select only genus and the existing b_columns
filtered_df = filtered_df[["genus"] + existing_b_columns]

# Group by genus (no aggregation of b_columns)
grouped = filtered_df.groupby("genus")

# Example: count rows per genus (replace with other agg if needed)
result = grouped.size().reset_index(name="count")

# Save to CSV
result.to_csv(output_file, index=False)

print(f"\nGrouped genus counts saved to:\n{output_file}")
