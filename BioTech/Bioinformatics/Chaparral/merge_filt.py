import os

def merge_filtered_fastq_files(base_path):
    # List barcode directories
    barcode_dirs = [os.path.join(base_path, d) for d in os.listdir(base_path)
                    if os.path.isdir(os.path.join(base_path, d)) and "barcode" in d.lower()]

    if not barcode_dirs:
        print("No barcode directories found.")
        return

    for barcode_dir in barcode_dirs:
        print(f"Merging filtered FASTQ files in {barcode_dir}...")

        # Find filtered fastq files in barcode directory
        filt_files = [os.path.join(barcode_dir, f) for f in os.listdir(barcode_dir)
                      if f.endswith("_trimmed_R_filt.fastq")]

        if not filt_files:
            print(f"  No filtered FASTQ files found in {barcode_dir}, skipping.")
            continue

        combined_path = os.path.join(barcode_dir, "combined_trimmed_R_filt.fastq")

        with open(combined_path, "wb") as outfile:
            for filt_file in filt_files:
                print(f"  Adding {os.path.basename(filt_file)}")
                with open(filt_file, "rb") as infile:
                    # Copy file content in binary mode
                    outfile.write(infile.read())

        print(f"  Combined FASTQ saved to {combined_path}\n")

if __name__ == "__main__":
    base_path = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\Round2"
    merge_filtered_fastq_files(base_path)
