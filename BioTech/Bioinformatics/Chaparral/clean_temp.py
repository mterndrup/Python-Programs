import os
import shutil

# Parent directory containing barcode folders
parent_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Wild2"

for foldername in os.listdir(parent_dir):
    barcode_path = os.path.join(parent_dir, foldername)
    if not os.path.isdir(barcode_path):
        continue

    # Delete all fasta files in the barcode folder
    for file in os.listdir(barcode_path):
        if file.endswith(".fasta"):
            fasta_path = os.path.join(barcode_path, file)
            try:
                os.remove(fasta_path)
                print(f"Deleted fasta file: {fasta_path}")
            except Exception as e:
                print(f"Error deleting {fasta_path}: {e}")

    # Delete ITSx folder if it exists in the barcode folder
    itsx_path = os.path.join(barcode_path, "ITSx")
    if os.path.isdir(itsx_path):
        try:
            shutil.rmtree(itsx_path)
            print(f"Deleted ITSx folder: {itsx_path}")
        except Exception as e:
            print(f"Error deleting ITSx folder {itsx_path}: {e}")

    # Go to Raw subfolder and delete cutadapt fastq files
    raw_path = os.path.join(barcode_path, "Raw")
    if os.path.isdir(raw_path):
        for file in os.listdir(raw_path):
            if file.endswith(".fastq") and "cutadapt" in file:
                fastq_path = os.path.join(raw_path, file)
                try:
                    os.remove(fastq_path)
                    print(f"Deleted fastq file: {fastq_path}")
                except Exception as e:
                    print(f"Error deleting {fastq_path}: {e}")
    else:
        print(f"No Raw folder in {barcode_path}, skipping Raw file deletion.")
