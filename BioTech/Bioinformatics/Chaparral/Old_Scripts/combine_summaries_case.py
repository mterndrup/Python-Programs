import pandas as pd


def update_csv_with_summary(summary_csv_path, csv_path, output_csv_path, barcode_col):
    # Step 1: Read summary counts from CSV file
    summary_df = pd.read_csv(summary_csv_path, encoding='utf-8')

    summary_counts = {}
    for _, row in summary_df.iterrows():
        species_str = str(row[0])  # First column: species/genus string
        try:
            count = int(row[1])  # Second column: count
        except (ValueError, TypeError):
            continue

        parts = species_str.strip().split()
        genus = parts[0].lower() if len(parts) > 0 else ''
        species = ' '.join(parts[1:]).lower() if len(parts) > 1 else ''
        key = (genus, species)

        # ✅ Accumulate counts instead of overwriting
        if key in summary_counts:
            summary_counts[key] += count
        else:
            summary_counts[key] = count

    print("📋 Parsed summary counts:")
    for k, v in summary_counts.items():
        print(f"  - {k}: {v}")

    # Step 2: Load master CSV and clean column names
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

    # Step 3: Ensure barcode-specific column exists
    if barcode_col not in df.columns:
        df[barcode_col] = 0
    else:
        df[barcode_col] = pd.to_numeric(df[barcode_col], errors='coerce').fillna(0).astype(int)

    # Step 4: Update existing or append new rows
    for (genus, species), count in summary_counts.items():
        print(f"🔍 Looking for: Genus='{genus}', Species='{species}'")
        match = (df['Genus'] == genus) & (df['Species'] == species)
        if match.any():
            print(f"✅ Match found for {genus} {species}: adding {count}")
            df.loc[match, barcode_col] += count
        else:
            print(f"➕ No match found for {genus} {species}, adding new row")
            new_row = {
                'Genus': genus,
                'Species': species,
                'Count': 0,
                barcode_col: count
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


# === Run the script ===
if __name__ == "__main__":
    barcode = input("Enter the barcode number (e.g., 21): ").strip()
    barcode_col = f"b{barcode}"

    base_path = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                 r"\Bioinformatics\Chaparral")

    summary_csv = (
        fr"{base_path}\wildfirePlants-DNA-nanopore-sequence\fastq_pass\barcode{barcode}\barcode{barcode}_sum.csv"
    )
    csv_file = fr"{base_path}\dna_counts.csv"
    output_csv = fr"{base_path}\dna_counts_updated_barcode{barcode}.csv"

    update_csv_with_summary(summary_csv, csv_file, output_csv, barcode_col)
