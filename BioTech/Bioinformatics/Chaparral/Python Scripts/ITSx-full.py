import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# === CONFIGURATION ===

# Original FASTA path (WSL format)
input_fasta = "/mnt/c/Users/ketgl/OneDrive/Documents/Glowing-Fungi/assembly.fasta"

# Path to ITSx executable
itsx_path = "/mnt/c/Users/ketgl/Downloads/ITSx_1.1.3/ITSx_1.1.3/ITSx"

# Number of CPUs
cpus = "4"

# Maximum allowed contig length for ITSx (HMMER limit is 100K)
max_len = 99000

# Output paths
base_dir = os.path.dirname(input_fasta)
split_fasta = os.path.join(base_dir, "assembly_contigs_split.fasta")
itsx_output_dir = os.path.join(base_dir, "ITSx")
os.makedirs(itsx_output_dir, exist_ok=True)
output_prefix = os.path.join(itsx_output_dir, "assembly_ITSx")

# === STEP 1: SPLIT LONG CONTIGS ===

print("Splitting long contigs...")

split_records = []
for record in SeqIO.parse(input_fasta, "fasta"):
    seq = record.seq
    if len(seq) <= max_len:
        # Keep short contigs as-is
        split_records.append(record)
    else:
        # Split long contigs into chunks
        for i in range(0, len(seq), max_len):
            chunk = seq[i:i + max_len]
            chunk_id = f"{record.id}_part_{i // max_len + 1}"
            split_records.append(SeqRecord(chunk, id=chunk_id, description=""))

# Write new FASTA
SeqIO.write(split_records, split_fasta, "fasta")
print(f"Contigs written to: {split_fasta}")

# === STEP 2: RUN ITSx ===

cmd = [
    "perl",
    itsx_path,
    "-i", split_fasta,
    "-o", output_prefix,
    "--cpu", cpus,
    "--preserve", "T",
    "--allow_partial",
    "--save_regions", "all"
]

print("Running ITSx...")
try:
    subprocess.run(cmd, check=True)
    print("ITSx completed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error running ITSx: {e}")
