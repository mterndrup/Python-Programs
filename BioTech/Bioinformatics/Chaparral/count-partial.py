import os
import sys
import csv
import glob

# Get barcode input
barcode = input("Enter the barcode number (e.g., 05 for barcode05): ").strip().zfill(2)
primertype = "FITS"

# Build directory path
input_dir = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
             r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
             fr"\fastq_pass\barcode{barcode}")

folder_name = os.path.basename(input_dir.rstrip(os.sep))

# Find latest blast log file
log_files = sorted(
    glob.glob(os.path.join(input_dir, f"{folder_name}_{primertype}_blast_log_*.txt")),
    key=os.path.getmtime
)

if not log_files:
    sys.exit(f"❌ No BLAST log files found in: {input_dir}")

log_file_path = log_files[-1]
print(f"📄 Using log file: {log_file_path}")

summary_csv_path = os.path.join(input_dir, f"{folder_name}_sum.csv")
species_counts = {}

# Read and process the log file
with open(log_file_path, "r", encoding="utf-8") as log_file:
    header = log_file.readline()  # skip header
    for line in log_file:
        parts = line.strip().split("\t")
        if len(parts) < 4:
            continue
        _, status, _, top_hit = parts
        if status != "Included":
            continue
        if "~" in top_hit:
            species = top_hit.split("~")[0].strip()
            species_counts[species] = species_counts.get(species, 0) + 1

# Write summary CSV
with open(summary_csv_path, "w", newline='', encoding="utf-8") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=["Species", "Count"])
    writer.writeheader()
    for species, count in sorted(species_counts.items(), key=lambda x: x[1], reverse=True):
        writer.writerow({"Species": species, "Count": count})

print(f"\n📊 Summary CSV saved: {summary_csv_path}")
print("✅ Summary generated from log file.")
