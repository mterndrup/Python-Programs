import pandas as pd

def update_csv_with_summary(summary_txt_path, csv_path, output_csv_path):
    # Step 1: Read summary counts from tab-separated text file
    summary_counts = {}
    with open(summary_txt_path, 'r', encoding='utf-8') as f:
        header_skipped = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if not header_skipped:
                header_skipped = True
                continue  # skip the header
            try:
                species_str, count_str = line.split('\t')
                count = int(count_str.strip())
            except ValueError:
                continue

            parts = species_str.strip().split()
            if len(parts) < 2:
                genus = parts[0].strip().lower()
                species = ''
            else:
                genus = parts[0].strip().lower()
                species = ' '.join(parts[1:]).strip().lower()
            summary_counts[(genus, species)] = count

    print("📋 Parsed summary counts:")
    for k, v in summary_counts.items():
        print(f"  - {k}: {v}")

    # Step 2: Load CSV and clean column names
    df = pd.read_csv(csv_path, sep=',', encoding='utf-8')
    df.columns = df.columns.str.strip().str.replace('\ufeff', '')

    print("\n🔍 Detected CSV columns:")
    for col in df.columns:
        print(f"  - '{col}'")

    if 'Genus' not in df.columns or 'Species' not in df.columns:
        raise ValueError(
            f"CSV must contain 'Genus' and 'Species' columns.\nFound columns: {list(df.columns)}"
        )

    # Normalize existing Genus/Species in CSV
    df['Genus'] = df['Genus'].astype(str).str.strip().str.lower()
    df['Species'] = df['Species'].astype(str).str.strip().str.lower()

    # Step 3: Ensure b21 column exists
    if 'b21' not in df.columns:
        df['b21'] = 0
    else:
        df['b21'] = pd.to_numeric(df['b21'], errors='coerce').fillna(0).astype(int)

    # Step 4: Update existing or append new rows
    for (genus, species), count in summary_counts.items():
        print(f"🔍 Looking for: Genus='{genus}', Species='{species}'")
        match = (df['Genus'] == genus) & (df['Species'] == species)
        if match.any():
            print(f"✅ Match found for {genus} {species}: adding {count}")
            df.loc[match, 'b21'] += count
        else:
            print(f"➕ No match found for {genus} {species}, adding new row")
            new_row = {
                'Genus': genus,
                'Species': species,
                'Count': 0,
                'b21': count
            }
            for col in df.columns:
                if col not in new_row:
                    new_row[col] = '' if df[col].dtype == object else 0
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

    # Step 5: Sort by Genus then Species
    df = df.sort_values(by=['Genus', 'Species'], ascending=[True, True], ignore_index=True)

    # Step 6: Save to new file
    df.to_csv(output_csv_path, sep=',', index=False, encoding='utf-8')
    print(f"\n✅ Updated CSV saved to:\n{output_csv_path}")

# === Run it ===
if __name__ == "__main__":
    summary_txt = (
        r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
        r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
        r"\fastq_pass\barcode21\summary_counts.txt"
    )
    csv_file = (
        r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
        r"\Bioinformatics\Chaparral\dna_counts.csv"
    )
    output_csv = (
        r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
        r"\Bioinformatics\Chaparral\dna_counts_updated.csv"
    )

    update_csv_with_summary(summary_txt, csv_file, output_csv)
