import os
import csv

# Paths
base_dir = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\Round2"
logs_dir = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\Logs"

# ITSx region types
regions = ['5_8S', 'ITS1', 'ITS2', 'LSU', 'SSU', 'full', 'no_detections']

def count_fasta_sequences(fasta_path):
    count = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count

def convert_barcode_to_prefix(barcode_folder):
    """Convert barcode folder name like 'barcode01' to 'b01'"""
    number = barcode_folder.replace("barcode", "")
    return f"b{number}"

def main():
    os.makedirs(logs_dir, exist_ok=True)

    csv_file = os.path.join(logs_dir, 'ITSx_counts_summary.csv')
    log_file = os.path.join(logs_dir, 'ITSx_counts_log.txt')

    with open(csv_file, 'w', newline='') as csvfile, open(log_file, 'w') as log:
        fieldnames = ['Sample'] + regions
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for barcode_folder in os.listdir(base_dir):
            barcode_path = os.path.join(base_dir, barcode_folder)
            if os.path.isdir(barcode_path) and barcode_folder.startswith("barcode"):

                sample_prefix = convert_barcode_to_prefix(barcode_folder)
                itsx_dir = os.path.join(barcode_path, "ITSx")
                counts = {region: 0 for region in regions}

                for region in regions:
                    if region == 'no_detections':
                        filename = f"{sample_prefix}_combined_trimmed_R_filt_ITSx_no_detections.fasta"
                    else:
                        filename = f"{sample_prefix}_combined_trimmed_R_filt_ITSx.{region}.fasta"

                    fasta_path = os.path.join(itsx_dir, filename)

                    if os.path.isfile(fasta_path):
                        try:
                            counts[region] = count_fasta_sequences(fasta_path)
                        except Exception as e:
                            log.write(f"Error reading {fasta_path}: {e}\n")
                    else:
                        log.write(f"File not found: {fasta_path}\n")

                writer.writerow({'Sample': barcode_folder, **counts})
                log.write(f"{barcode_folder}: {counts}\n")

    print(f"✅ Summary saved: {csv_file}")
    print(f"📝 Log written: {log_file}")

if __name__ == "__main__":
    main()
