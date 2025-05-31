import subprocess
import os

# === Set vsearch path ===
vsearch_path = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                r"\Bioinformatics\Chaparral\vsearch-2.30.0-win-x86_64"
                r"\vsearch-2.30.0-win-x86_64\bin\vsearch.exe")

# === Input paths ===
input_its2 = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
              r"\Bioinformatics\Tribal-Medicine\Cacti\DNA_Sequences\barcode33"
              r"\barcode33_PITS.fasta")

# === Setup output paths ===
base_dir = os.path.dirname(input_its2)
base_name = os.path.splitext(os.path.basename(input_its2))[0]

derep_fasta = os.path.join(base_dir, f"{base_name}_derep.fasta")
centroids_fasta = os.path.join(base_dir, f"{base_name}_centroids.fasta")
uc_file = os.path.join(base_dir, f"{base_name}_clusters.uc")
sorted_fasta = os.path.join(base_dir, f"{base_name}_sorted_centroids.fasta")
dominant_fasta = os.path.join(base_dir, f"{base_name}_dominant.fasta")
mapping_output = os.path.join(base_dir, "dominant_vs_full.tsv")

# === Step 1: Dereplicate ===
print("Step 1: Dereplicating ITS2 sequences...")
subprocess.run([
    vsearch_path,
    "--derep_fulllength", input_its2,
    "--output", derep_fasta,
    "--sizeout"
], check=True)

# === Step 2: Cluster at 99% ===
print("Step 2: Clustering sequences at 99% identity...")
subprocess.run([
    vsearch_path,
    "--cluster_fast", derep_fasta,
    "--id", "0.99",
    "--strand", "both",
    "--centroids", centroids_fasta,
    "--uc", uc_file,
    "--sizeout"
], check=True)


# === Step 4: Extract top (dominant) sequence ===
print("Step 4: Extracting the most abundant sequence...")
top_seq = ""
top_seq_id = ""
with open(centroids_fasta, "r") as f:
    for line in f:
        if line.startswith(">"):
            top_seq_id = line.strip()
            continue
        top_seq = line.strip()
        break

with open(dominant_fasta, "w") as out:
    out.write(f"{top_seq_id}\n{top_seq}\n")
