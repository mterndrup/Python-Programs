import os
import csv

# Paths
base_dir = r"/BioTech/Bioinformatics/Chaparral/Round2"
regions = ['5_8S', 'ITS1', 'ITS2', 'LSU', 'SSU', 'full', 'no_detections']

def convert_barcode_to_prefix(barcode_folder):
    number = barcode_folder.replace("barcode", "")
    return f"b{number}"

def parse_fasta_headers(fasta_path):
    """Returns a set of sequence IDs from the given FASTA file."""
    ids = set()
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]  # Keep only the first word (ID)
                ids.add(seq_id)
    return ids

def main():
    for barcode_folder in os.listdir(base_dir):
        barcode_path = os.path.join(base_dir, barcode_folder)
        if os.path.isdir(barcode_path) and barcode_folder.startswith("barcode"):

            print(f"🔍 Processing {barcode_folder}...")

            sample_prefix = convert_barcode_to_prefix(barcode_folder)
            itsx_dir = os.path.join(barcode_path, "ITSx-v1")

            seq_region_map = {}  # {seqID: {region: 0/1}}

            for region in regions:
                if region == 'no_detections':
                    filename = f"{sample_prefix}_combined_trimmed_R_filt_ITSx_no_detections.fasta"
                else:
                    filename = f"{sample_prefix}_combined_trimmed_R_filt_ITSx.{region}.fasta"

                fasta_path = os.path.join(itsx_dir, filename)

                if os.path.isfile(fasta_path):
                    try:
                        ids = parse_fasta_headers(fasta_path)
                        for seq_id in ids:
                            if seq_id not in seq_region_map:
                                seq_region_map[seq_id] = {r: 0 for r in regions}
                            seq_region_map[seq_id][region] = 1
                    except Exception as e:
                        print(f"[ERROR] Couldn't read {fasta_path}: {e}")
                else:
                    print(f"[WARN] Missing: {fasta_path}")

            # Write output CSV to the barcode folder
            output_csv = os.path.join(barcode_path, f"{barcode_folder}_ITSx_seq_region_counts.csv")
            with open(output_csv, 'w', newline='') as csvfile:
                fieldnames = ['seqID'] + regions
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()

                for seq_id, presence in seq_region_map.items():
                    row = {'seqID': seq_id, **presence}
                    writer.writerow(row)

            print(f"✅ Region count CSV saved: {output_csv}\n")

if __name__ == "__main__":
    main()
