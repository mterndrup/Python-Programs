import pandas as pd
import os

genus_type_map = {
    "Mucor": "Heat-Activated Germinators",
    "Thermomyces": "Heat-Activated Germinators",
    "Rhizopus": "Heat-Activated Germinators",
    "Penicillium": "Rapid Colonizers",
    "Aspergillus": "Rapid Colonizers",
    "Fusarium": "Rapid Colonizers",
    "Trichoderma": "Rapid Colonizers",
    "Daldinia": "PyOM Degraders",
    "Cryomyces": "PyOM Degraders",
    "Lenzites": "PyOM Degraders",
    "Ganoderma": "Plant Debris Decomposers",
    "Colletotrichum": "Plant Debris Decomposers",
    "Didymella": "Plant Debris Decomposers",
    "Henningsomyces": "Plant Debris Decomposers",
    "Peziza": "Plant Debris Decomposers",
    "Boletus": "Plant Debris Decomposers",
    "Pisolithus": "Mycorrhizal Establishers",
    "Serendipita": "Mycorrhizal Establishers",
    "Amanita": "Mycorrhizal Establishers",
}


def filter_and_average_count(file_path, output_path):
    # Read CSV (default delimiter is comma)
    df = pd.read_csv(file_path)

    # Strip whitespace from column headers
    df.columns = df.columns.str.strip()

    print("Columns found in CSV:", df.columns.tolist())  # debug line — remove if you want

    # Confirm required columns exist
    if 'Genus' not in df.columns or 'Species' not in df.columns or 'Count' not in df.columns:
        raise ValueError("CSV file must contain 'Genus', 'Species', and 'Count' columns.")

    # Normalize Genus and Species columns
    df['Genus'] = df['Genus'].str.capitalize()
    df['Species'] = df['Species'].str.lower()

    # Filter by genus list
    filtered_df = df[df['Genus'].isin(genus_type_map.keys())].copy()

    # Exclude rows with Species == 'sp.' or contains 'uncultured'
    filtered_df = filtered_df[~filtered_df['Species'].str.contains(r'(^sp\.?$|uncultured)', regex=True)]

    # Group by Genus and average only the Count column
    averaged_df = filtered_df.groupby('Genus')['Count'].mean().reset_index()

    # Map Type to each Genus
    averaged_df['Type'] = averaged_df['Genus'].map(genus_type_map)

    # Reorder columns to Genus, Type, Count
    averaged_df = averaged_df[['Genus', 'Type', 'Count']]

    # Save to CSV
    averaged_df.to_csv(output_path, index=False)
    print(f"Averaged Count data saved to: {output_path}")


if __name__ == "__main__":
    base_dir = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral"
    input_file = os.path.join(base_dir, "dna_counts.csv")
    output_file = os.path.join(base_dir, "averaged_count_dna_counts.csv")

    filter_and_average_count(input_file, output_file)
