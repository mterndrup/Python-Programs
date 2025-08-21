import os
import sys
import datetime

# Base directory
BASE_DIR = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox"
LOG_DIR = os.path.join(BASE_DIR, "Logs")

# Ensure log directory exists
os.makedirs(LOG_DIR, exist_ok=True)

# Script info
SCRIPT_NAME = os.path.basename(__file__)

# Timestamped log file
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = os.path.join(LOG_DIR, f"log_{timestamp}.txt")

def log(message):
    """Write to terminal and log file."""
    with open(log_file, "a") as f:
        f.write(message + "\n")
    print(message)

def count_fasta_sequences(file_path):
    """Count sequences in a FASTA file by counting '>' lines."""
    count = 0
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count

def extract_primer_type(filename):
    """Extract primer type (18S, FITS, PITS) from filename."""
    fname = filename.lower()
    if "18s" in fname:
        return "18S"
    elif "fits" in fname:
        return "FITS"
    elif "pits" in fname:
        return "PITS"
    else:
        return "OTHER"

def main():
    start_time = datetime.datetime.now()
    log("="*60)
    log(f"Script: {SCRIPT_NAME}")
    log(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log("="*60)

    total_sequences_all = 0
    barcode_totals = {}
    primer_grand_totals = {"18S": 0, "FITS": 0, "PITS": 0, "OTHER": 0}

    # Only scan Cactus and Cactus1
    for folder in ["Cactus", "Cactus1"]:
        folder_path = os.path.join(BASE_DIR, folder)
        if not os.path.exists(folder_path):
            log(f"Skipping missing folder: {folder_path}")
            continue

        log(f"\nProcessing folder: {folder_path}")

        for barcode in os.listdir(folder_path):
            barcode_path = os.path.join(folder_path, barcode)
            if not os.path.isdir(barcode_path):
                continue

            log(f"\n  Barcode: {barcode}")
            barcode_total = 0

            for primer_file in os.listdir(barcode_path):
                primer_path = os.path.join(barcode_path, primer_file)

                # ✅ Only process FASTA-like files (allow extensionless)
                if not os.path.isfile(primer_path):
                    continue
                if not (primer_file.lower().endswith((".fa", ".fasta")) or "." not in primer_file):
                    continue

                seq_count = count_fasta_sequences(primer_path)
                barcode_total += seq_count

                # ✅ add to primer grand total (collapsed by type)
                primer_type = extract_primer_type(primer_file)
                primer_grand_totals[primer_type] += seq_count

                log(f"    Primer: {primer_file} -> {seq_count} sequences")

            barcode_totals[barcode] = barcode_total
            total_sequences_all += barcode_total

            log(f"  Total sequences for {barcode}: {barcode_total}")

    log("\n" + "="*60)
    log("SUMMARY")
    for barcode, total in barcode_totals.items():
        log(f"  {barcode}: {total} sequences")

    log(f"\nGrand total sequences (all barcodes): {total_sequences_all}")
    log(f"Grand total 18S: {primer_grand_totals['18S']}")
    log(f"Grand total FITS: {primer_grand_totals['FITS']}")
    log(f"Grand total PITS: {primer_grand_totals['PITS']}")
    if primer_grand_totals["OTHER"] > 0:
        log(f"Grand total OTHER: {primer_grand_totals['OTHER']}")

    log("="*60)

    end_time = datetime.datetime.now()
    duration = end_time - start_time
    log(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"Total processing time: {duration}")
    log("="*60)

if __name__ == "__main__":
    main()
