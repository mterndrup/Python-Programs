import os
import shutil

# Paths
source_root = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Wild2"
dest_root = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Wild4"

# Ensure destination exists
os.makedirs(dest_root, exist_ok=True)

# Loop through all folders in source_root
for item in os.listdir(source_root):
    source_item_path = os.path.join(source_root, item)

    # Only process directories that start with "barcode"
    if os.path.isdir(source_item_path) and item.lower().startswith("barcode"):
        dest_item_path = os.path.join(dest_root, item)
        os.makedirs(dest_item_path, exist_ok=True)

        # Copy Raw folder if it exists
        raw_source = os.path.join(source_item_path, "Raw")
        raw_dest = os.path.join(dest_item_path, "Raw")
        if os.path.isdir(raw_source):
            shutil.copytree(raw_source, raw_dest, dirs_exist_ok=True)
            print(f"Copied Raw folder for {item}")

        # Copy specific fasta files
        fasta_suffixes = ["_18S.fasta", "_FITS.fasta", "_PITS.fasta"]
        for suffix in fasta_suffixes:
            fasta_source = os.path.join(source_item_path, f"{item}{suffix}")
            if os.path.isfile(fasta_source):
                shutil.copy2(fasta_source, dest_item_path)
                print(f"Copied {os.path.basename(fasta_source)} for {item}")

print("Copy complete.")
