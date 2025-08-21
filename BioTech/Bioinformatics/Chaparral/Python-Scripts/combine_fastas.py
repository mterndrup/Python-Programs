import os
from Bio import SeqIO
from datetime import datetime
import sys

# Parent directory where barcode folders live
parent_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus2"

# Log directory
log_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Logs"
os.makedirs(log_dir, exist_ok=True)

# Timestamped log file
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_path = os.path.join(log_dir, f"combine_fasta_run_{timestamp}.log")

# Get the current script name
script_name = os.path.basename(sys.argv[0])

def log_print(message, log_file):
    print(message)
    log_file.write(message + "\n")

start_time = datetime.now()

with open(log_path, 'a') as log_file:
    log_print(f"=== {script_name} Run Started: {start_time} ===", log_file)

    total_barcodes = 0
    for foldername in sorted(os.listdir(parent_dir)):
        if foldername.startswith("barcode"):
            barcode_folder = os.path.join(parent_dir, foldername)
            if not os.path.isdir(barcode_folder):
                continue

            # Only select fasta files ending with _18S.fasta, _FITS.fasta, _PITS.fasta
            fasta_files = [
                f for f in os.listdir(barcode_folder)
                if f.endswith("_18S.fasta") or f.endswith("_FITS.fasta") or f.endswith("_PITS.fasta")
            ]

            if not fasta_files:
                log_print(f"No 18S/FITS/PITS fasta files found in {foldername}, skipping.", log_file)
                continue

            combined_fasta_path = os.path.join(barcode_folder, f"{foldername}_combined.fasta")
            seq_ids_seen = set()
            count_written = 0

            with open(combined_fasta_path, 'w') as out_handle:
                for fasta_file in fasta_files:
                    fasta_path = os.path.join(barcode_folder, fasta_file)
                    log_print(f"Processing {fasta_path}", log_file)

                    for record in SeqIO.parse(fasta_path, "fasta"):
                        original_id = record.id
                        new_id = original_id
                        suffix = 1
                        while new_id in seq_ids_seen:
                            new_id = f"{original_id}_{suffix}"
                            suffix += 1
                        record.id = new_id
                        record.description = ""
                        seq_ids_seen.add(new_id)
                        SeqIO.write(record, out_handle, "fasta")
                        count_written += 1

            log_print(f"{foldername}: {count_written} sequences written to {combined_fasta_path}", log_file)
            total_barcodes += 1

    end_time = datetime.now()
    elapsed = end_time - start_time

    log_print(f"\n=== {script_name} Run Completed: {end_time} ===", log_file)
    log_print(f"Total processing time: {elapsed}", log_file)
    log_print(f"Total barcodes processed: {total_barcodes}", log_file)

print(f"Combined fasta files are created in each barcode folder.")
print(f"Log saved to: {log_path}")
