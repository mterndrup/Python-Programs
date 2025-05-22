import os
import gzip
import shutil

base_dir = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence\fastq_pass"

for barcode_folder in sorted(os.listdir(base_dir)):
    barcode_path = os.path.join(base_dir, barcode_folder)

    if not os.path.isdir(barcode_path) or not barcode_folder.startswith("barcode"):
        continue

    print(f"📂 Processing {barcode_folder}...")

    for item in os.listdir(barcode_path):
        if item.endswith(".gz"):
            gz_path = os.path.join(barcode_path, item)
            out_file = os.path.splitext(gz_path)[0]  # remove .gz

            if os.path.exists(out_file):
                print(f"  ⚠️ Skipping {item}, {os.path.basename(out_file)} already exists.")
                continue

            print(f"  🗜 Decompressing {item}...")
            with gzip.open(gz_path, 'rb') as f_in:
                with open(out_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

print("✅ All .gz files decompressed.")