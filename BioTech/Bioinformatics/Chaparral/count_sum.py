import pandas as pd

def fill_and_sort_counts(csv_path, output_path):
    # Load CSV
    df = pd.read_csv(csv_path, sep=',', encoding='utf-8')
    df.columns = df.columns.str.strip().str.replace('\ufeff', '')

    # Detect b-columns
    b_columns = [col for col in df.columns if col.startswith('b')]
    if not b_columns:
        raise ValueError("No columns starting with 'b' were found.")
    print(f"✅ Detected b-columns: {b_columns}")

    # Replace NaNs and convert to numeric
    df[b_columns] = df[b_columns].fillna(0)
    df[b_columns] = df[b_columns].apply(pd.to_numeric, errors='coerce').fillna(0)

    # Fill 'Count' with the number of b-columns with values > 0
    df['Count'] = (df[b_columns] > 0).sum(axis=1)

    # Sort: First by Count (descending), then Genus (A–Z), then Species (A–Z)
    df.sort_values(by=['Count', 'Genus', 'Species'], ascending=[False, True, True], inplace=True)

    # Save updated CSV
    df.to_csv(output_path, index=False, sep=',', encoding='utf-8')
    print(f"✅ Sorted CSV saved to: {output_path}")


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

    fill_and_sort_counts(input_csv, output_csv)
