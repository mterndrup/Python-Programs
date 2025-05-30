import os
import subprocess

# Path to NanoFilt executable (adjust as needed)
nanofilt_path = r"C:\Users\ketgl\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0\LocalCache\local-packages\Python312\Scripts\NanoFilt.exe"

# Base directory containing barcode folders
base_dir = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
            r"\Bioinformatics\Tribal-Medicine\Cacti\DNA_Sequences")

# Get all barcode directories (folders starting with "barcode")
barcode_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and d.startswith("barcode")]

for barcode in barcode_dirs:
    raw_folder = os.path.join(base_dir, barcode, "Raw")
    if not os.path.exists(raw_folder):
        print(f"⚠️ Raw folder not found for {barcode}, skipping.")
        continue

    print(f"\n🔍 Processing {barcode}...")

    # Find all trimmed fastq files in the Raw folder
    fastq_files = [f for f in os.listdir(raw_folder) if f.endswith("_trimmed.fastq")]

    for fq_file in fastq_files:
        fq_path = os.path.join(raw_folder, fq_file)
        output_path = fq_path.replace(".fastq", "_filtered.fastq")

        print(f"  🧪 Filtering {fq_file} -> {output_path}")

        try:
            with open(fq_path, "r") as in_f, open(output_path, "w") as out_f:
                subprocess.run([nanofilt_path, "-q", "10", "-l", "200"], stdin=in_f, stdout=out_f, check=True)
        except subprocess.CalledProcessError as e:
            print(f"  ❌ Error filtering {fq_file}: {e}")
        except FileNotFoundError:
            print(f"  ❌ NanoFilt not found at {nanofilt_path}. Please verify the path.")

    # Delete any zero-byte files in the Raw folder after processing this barcode folder
    for filename in os.listdir(raw_folder):
        file_path = os.path.join(raw_folder, filename)
        if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
            os.remove(file_path)
            print(f"  🗑️ Deleted zero-byte file: {filename}")

print("\n✅ Done filtering all _trimmed.fastq files.")
