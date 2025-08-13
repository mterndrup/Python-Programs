import os
import pandas as pd

# Define the path to the directory containing barcode folders
base_dir = "../Wildfire-Sequences/fastq_pass"
barcodes = [f"barcode{str(i).zfill(2)}" for i in range(1, 65)]

# Define the required columns and fallbacks
required_columns = ["genus", "species", "vsearch_count", "dada2_count", "only_in"]
fallback_column_map = {
    "genus": ["genus", "Genus", "GENUS"],
    "species": ["species", "Species", "SPECIES"],
    "vsearch_count": ["vsearch_count", "vsearch", "VSEARCH_COUNT", "vsearchcount"],
    "dada2_count": ["dada2_count", "dada2", "DADA2_COUNT", "dada2count"],
    "only_in": ["only_in", "Only_in", "ONLY_IN"]
}

# Store all rows for later merging
all_data = []

for barcode in barcodes:
    filename = f"comparison_vsearch_dada2_{barcode[7:]}.csv"
    filepath = os.path.join(base_dir, barcode, filename)

    if not os.path.isfile(filepath):
        print(f"CSV not found in {barcode}: {filename}")
        continue

    print(f"Reading: {filepath}")
    try:
        df = pd.read_csv(filepath)
        df.columns = df.columns.str.strip().str.lower()

        # Map columns using the fallback map
        mapped_columns = {}
        for col_key, aliases in fallback_column_map.items():
            for alias in aliases:
                alias_clean = alias.strip().lower()
                if alias_clean in df.columns:
                    mapped_columns[col_key] = alias_clean
                    break

        if not all(key in mapped_columns for key in required_columns):
            print(f"Skipping {filepath}: missing required columns. Found: {list(df.columns)}")
            continue

        # Rename to standard column names
        df = df.rename(columns={v: k for k, v in mapped_columns.items()})

        # Keep only necessary columns
        df = df[["genus", "species", "vsearch_count", "dada2_count"]]
        df["barcode"] = barcode

        all_data.append(df)

    except Exception as e:
        print(f"Error processing {filepath}: {e}")

# Merge all data
if all_data:
    combined = pd.concat(all_data, ignore_index=True)

    # Group by genus/species, sum counts, collect barcode list
    grouped = combined.groupby(["genus", "species"]).agg({
        "vsearch_count": "sum",
        "dada2_count": "sum",
        "barcode": lambda x: ",".join(sorted(set(x)))
    }).reset_index()

    # Sort by dada2_count and vsearch_count descending
    grouped = grouped.sort_values(by=["dada2_count", "vsearch_count"], ascending=False)

    # Save to CSV
    output_path = os.path.join(base_dir, "species_summary_by_counts.csv")
    grouped.to_csv(output_path, index=False)
    print(f"Summary saved to: {output_path}")
else:
    print("No matching or usable CSV files found.")
