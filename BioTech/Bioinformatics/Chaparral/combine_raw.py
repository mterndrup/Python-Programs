import os
from Bio import SeqIO

def get_suffix(filename):
    # Remove extension first
    name = os.path.splitext(filename)[0]
    # If file is .fastq.gz, remove .gz as well
    if filename.endswith(".fastq.gz"):
        name = os.path.splitext(name)[0]
    parts = name.split('_')
    # Assuming suffix always consists of last 4 parts
    suffix_parts = parts[-4:]
    suffix = "_" + "_".join(suffix_parts)
    return suffix

base_dir = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence\fastq_pass"

barcode_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and d.startswith("barcode")]

for barcode in barcode_dirs:
    raw_folder = os.path.join(base_dir, barcode, "Raw")
    if not os.path.exists(raw_folder):
        print(f"⚠️ Raw folder not found for {barcode}, skipping.")
        continue

    output_fasta = os.path.join(base_dir, barcode, f"{barcode}_combined.fasta")
    print(f"\nProcessing {barcode}...")
    print(f"Saving combined FASTA to {output_fasta}")

    # Only read files ending with '_trimmed_filtered.fastq'
    fastq_files = [f for f in os.listdir(raw_folder) if f.endswith("_trimmed_filtered.fastq")]

    with open(output_fasta, "w") as fasta_out:
        for fq_file in fastq_files:
            fq_path = os.path.join(raw_folder, fq_file)
            suffix = get_suffix(fq_file)
            print(f"  Reading {fq_file} with suffix '{suffix}'...")

            for record in SeqIO.parse(fq_path, "fastq"):
                header = f"@{record.id}{suffix}"
                fasta_out.write(f">{header}\n{str(record.seq)}\n")

print("\n✅ All done!")
