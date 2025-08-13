import os
import csv

# Base directory containing barcode folders
base_dir = r"/BioTech/Bioinformatics/Chaparral/Round2"

def read_fasta_ids(fasta_path):
    """Read sequence IDs from a FASTA file."""
    ids = set()
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line.split()[0].lstrip(">")
                ids.add(seq_id)
    return ids

def has_its_region(row):
    """Check if any ITS region (except no_detections) has a value > 0."""
    # Columns to check (all except 'seqID' and 'no_detections')
    its_columns = [col for col in row.keys() if col not in ('seqID', 'no_detections')]
    for col in its_columns:
        try:
            if int(row[col]) > 0:
                return True
        except ValueError:
            # Handle cases where conversion fails gracefully
            continue
    return False

def process_barcode_folder(barcode_folder):
    barcode_path = os.path.join(base_dir, barcode_folder)
    sample_prefix = barcode_folder.replace("barcode", "b")

    fasta_file = os.path.join(barcode_path, f"{sample_prefix}_combined_trimmed_R_filt_kingdom.fasta")
    csv_file = os.path.join(barcode_path, f"{barcode_folder}_ITSx_seq_region_counts.csv")
    output_csv = os.path.join(barcode_path, f"{barcode_folder}_ITSx_seq_region_filtered.csv")

    if not os.path.isfile(fasta_file):
        print(f"⚠️ Missing FASTA file: {fasta_file}")
        return
    if not os.path.isfile(csv_file):
        print(f"⚠️ Missing CSV file: {csv_file}")
        return

    fasta_ids = read_fasta_ids(fasta_file)

    with open(csv_file, 'r', newline='') as infile:
        reader = csv.DictReader(infile)  # comma-separated by default
        if 'seqID' not in reader.fieldnames:
            print(f"⚠️ Invalid or missing 'seqID' in CSV header: {csv_file}")
            return

        fieldnames = reader.fieldnames
        rows_to_write = []

        for row in reader:
            if row['seqID'] in fasta_ids and has_its_region(row):
                rows_to_write.append(row)

    if not rows_to_write:
        print(f"⚠️ No ITS regions found in: {csv_file}")
        return

    with open(output_csv, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)  # comma-separated default
        writer.writeheader()
        writer.writerows(rows_to_write)

    print(f"✅ Filtered CSV created: {output_csv}")

def main():
    for barcode_folder in sorted(os.listdir(base_dir)):
        barcode_path = os.path.join(base_dir, barcode_folder)
        if os.path.isdir(barcode_path) and barcode_folder.startswith("barcode"):
            process_barcode_folder(barcode_folder)

if __name__ == "__main__":
    main()
