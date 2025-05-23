import pandas as pd

# Mapping of Genus to Functional Group
GENUS_GROUP_MAP = {
    'Mucor': 'Heat-Activated Germinators',
    'Thermomyces': 'Heat-Activated Germinators',
    'Rhizopus': 'Heat-Activated Germinators',
    'Penicillum': 'Rapid Colonizers',
    'Aspergillus': 'Rapid Colonizers',
    'Fusarium': 'Rapid Colonizers',
    'Trichoderma': 'Rapid Colonizers',
    'Daldinia': 'PyOM Degraders',
    'Cryomyces': 'PyOM Degraders',
    'Lenzites': 'PyOM Degraders',
    'Ganoderma': 'Plant Debris Decomposers',
    'Colletotrichum': 'Plant Debris Decomposers',
    'Didymella': 'Plant Debris Decomposers',
    'Henningsomyces': 'Plant Debris Decomposers',
    'Peziza': 'Plant Debris Decomposers',
    'Boletus': 'Plant Debris Decomposers',
    'Pisolithus': 'Mycorrhizal Establishers',
    'Serendipita': 'Mycorrhizal Establishers',
    'Amanita': 'Mycorrhizal Establishers'
}


def fill_sort_and_summarize(csv_path, output_path, summary_path):
    # Load CSV
    df = pd.read_csv(csv_path, sep=',', encoding='utf-8')
    df.columns = df.columns.str.strip().str.replace('\ufeff', '')

    # Detect b-columns
    b_columns = [col for col in df.columns if col.startswith('b')]
    if not b_columns:
        raise ValueError("No columns starting with 'b' were found.")
    print(f"✅ Detected b-columns: {b_columns}")

    # Clean and process b-columns
    df[b_columns] = df[b_columns].fillna(0)
    df[b_columns] = df[b_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
    df['Count'] = (df[b_columns] > 0).sum(axis=1)

    # Sort full dataset
    df.sort_values(by=['Count', 'Genus', 'Species'], ascending=[False, True, True], inplace=True)
    df.to_csv(output_path, index=False, sep=',', encoding='utf-8')
    print(f"✅ Sorted CSV saved to: {output_path}")

    # Ensure Species is string for filtering
    df['Species'] = df['Species'].astype(str)

    # === Build summary ===
    summary_rows = []
    for genus, group in GENUS_GROUP_MAP.items():
        filtered = df[
            (df['Genus'] == genus) &
            (~df['Species'].str.lower().str.contains('uncultured')) &
            (~df['Species'].str.lower().str.contains(r'\bsp\b\.?'))
        ]
        if not filtered.empty:
            avg_count = filtered['Count'].mean()
            summary_rows.append({
                'Genus': genus,
                'Functional Group': group,
                'Average Count': round(avg_count, 2)
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(summary_path, index=False, sep=',', encoding='utf-8')
    print(f"✅ Genus summary saved to: {summary_path}")


# === Run It ===
if __name__ == "__main__":
    input_csv = (
        r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
        r"\Bioinformatics\Chaparral\dna_counts.csv"
    )
    output_csv = (
        r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
        r"\Bioinformatics\Chaparral\dna_counts_final_sorted.csv"
    )
    summary_csv = (
        r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
        r"\Bioinformatics\Chaparral\genus_functional_groups.csv"
    )

    fill_sort_and_summarize(input_csv, output_csv, summary_csv)
