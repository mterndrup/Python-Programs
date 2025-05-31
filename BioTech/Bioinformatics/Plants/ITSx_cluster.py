import subprocess
import os

# === Set vsearch path ===
vsearch_path = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                r"\Bioinformatics\Chaparral\vsearch-2.30.0-win-x86_64"
                r"\vsearch-2.30.0-win-x86_64\bin\vsearch.exe")

# === Input paths ===
input_its2 = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Tribal-Medicine\Cacti\DNA_Sequences\barcode34\ITSx\barcode34_PITS_ITSx.ITS2.fasta"
full_its_ref = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Tribal-Medicine\Cacti\DNA_Sequences\barcode34\ITSx\barcode34_PITS_ITSx.full.fasta"

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
    "--cluster_size", derep_fasta,
    "--id", "0.99",
    "--centroids", centroids_fasta,
    "--uc", uc_file,
    "--sizeout"
], check=True)

# === Step 3: Sort by size ===
print("Step 3: Sorting clusters by size...")
subprocess.run([
    vsearch_path,
    "--sortbysize", centroids_fasta,
    "--output", sorted_fasta
], check=True)

# === Step 4: Extract top (dominant) sequence ===
print("Step 4: Extracting the most abundant sequence...")
top_seq = ""
top_seq_id = ""
with open(sorted_fasta, "r") as f:
    for line in f:
        if line.startswith(">"):
            top_seq_id = line.strip()
            continue
        top_seq = line.strip()
        break

with open(dominant_fasta, "w") as out:
    out.write(f"{top_seq_id}\n{top_seq}\n")

print(f"Dominant ITS2 sequence saved to:\n{dominant_fasta}")

# === Step 5: Map dominant ITS2 to full ITS reference ===
print("Step 5: Mapping dominant ITS2 to full-length ITS reference...")
subprocess.run([
    vsearch_path,
    "--usearch_global", dominant_fasta,
    "--db", full_its_ref,
    "--id", "0.95",
    "--strand", "plus",
    "--blast6out", mapping_output
], check=True)

print(f"Mapping complete. Results saved to:\n{mapping_output}")

# === Step 6: Evaluate mapping quality ===
print("Step 6: Evaluating mapping quality...")
match_found = False
with open(mapping_output, "r") as f:
    for line in f:
        cols = line.strip().split("\t")
        if len(cols) >= 3:
            try:
                identity = float(cols[2])
                aln_len = int(cols[3])
                if identity >= 99.0 and aln_len >= 100:
                    match_found = True
                    break
            except ValueError:
                continue

if match_found:
    print("✅ Dominant ITS2 matches full ITS — likely valid.")
else:
    print("⚠️ No good match found. Full ITS may be contaminated or from a different source.")
