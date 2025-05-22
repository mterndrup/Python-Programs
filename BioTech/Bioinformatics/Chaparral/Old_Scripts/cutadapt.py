import os
import subprocess
import time
import urllib.error
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

# === Config ===
input_dir = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
             r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
             r"\fastq_pass\barcode05")
primertype = "FITS"
fwd_primer = "GGAAGTAAAAGTCGTAACAAGG"
rev_primer = "CAAGAGATCCGTTGTTGAAAGTT"
cutadapt_path = r"C:\Users\ketgl\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0\LocalCache\local-packages\Python312\Scripts\cutadapt.exe"

identity_threshold = 75  # Adjusted from 80% to 75%

folder_name = os.path.basename(input_dir.rstrip(os.sep))
input_files = sorted(f for f in os.listdir(input_dir) if f.endswith("_trimmed.fastq"))

if not input_files:
    print("No _trimmed.fastq files found.")
    exit()

final_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_final_trimmed_raw.fastq")

with open(final_output_path, "w") as final_out:
    for filename in input_files:
        input_path = os.path.join(input_dir, filename)
        temp_output_path = os.path.join(input_dir, f"temp_{filename}_cutadapt.fastq")
        cmd = [
            cutadapt_path, "-g", fwd_primer, "-a", rev_primer,
            "-e", "0.1", "--minimum-length", "60", "--overlap", "5",
            "--discard-untrimmed", "-o", temp_output_path, input_path
        ]
        print(f"🔪 Running Cutadapt on {filename}")
        try:
            subprocess.run(cmd, check=True)
            with open(temp_output_path, "r") as temp_in:
                final_out.write(temp_in.read())
            os.remove(temp_output_path)
        except subprocess.CalledProcessError as e:
            print(f"❌ Error running cutadapt on {filename}: {e}")

print(f"🎉 Trimmed FASTQ saved: {final_output_path}")

# === Convert to FASTA ===
fasta_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_final_trimmed.fasta")
with open(final_output_path, "r") as fastq_in, open(fasta_output_path, "w") as fasta_out:
    line_num = 0
    for line in fastq_in:
        line = line.strip()
        if line_num % 4 == 0:
            fasta_out.write(f">{line[1:]}\n")
        elif line_num % 4 == 1:
            fasta_out.write(line + "\n")
        line_num += 1

print(f"📂 FASTA saved: {fasta_output_path}")

# === Helper Functions ===
def blast_with_retry(seq, max_retries=3):
    wait_time = 5
    for attempt in range(max_retries):
        try:
            return NCBIWWW.qblast("blastn", "nt", seq, hitlist_size=1)
        except (urllib.error.URLError, ConnectionResetError) as e:
            print(f"🔁 Retry {attempt + 1}: {e}")
            time.sleep(wait_time)
            wait_time *= 2
    return None

def extract_genus_species(description):
    parts = description.split()
    if len(parts) >= 2:
        return parts[0], f"{parts[0]} {parts[1]}"
    elif parts:
        return parts[0], parts[0]
    return "unknown", "unknown"

# === Filtering ===
fungal_keywords = [
    "Fungi", "Ascomycota", "Basidiomycota", "Zygomycota", "Chytridiomycota",
    "Glomeromycota", "Microsporidia", "fungus", "mold", "yeast", "mycota", "mycetes"
]

# === BLAST & Output ===
fasta_blast_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_fungi_only.fasta")
log_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_blast_log.txt")

species_counts = {}

with open(fasta_output_path, "r") as fasta_in, \
     open(fasta_blast_output_path, "w") as fasta_out, \
     open(log_output_path, "w") as log_file:

    log_file.write("ID\tStatus\tReason\tTop_Hit\n")  # Log header

    for record in SeqIO.parse(fasta_in, "fasta"):
        seq = record.seq
        print(f"\n🔬 BLASTing: {record.id}")
        result_handle = blast_with_retry(seq)

        if result_handle is None:
            reason = "No BLAST result"
            print(f"  ❌ {reason}")
            log_file.write(f"{record.id}\tExcluded\t{reason}\tN/A\n")
            continue

        try:
            blast_record = NCBIXML.read(result_handle)
        except Exception as e:
            reason = f"BLAST parse error: {e}"
            print(f"  ❌ {reason}")
            log_file.write(f"{record.id}\tExcluded\t{reason}\tN/A\n")
            continue

        if blast_record.alignments and blast_record.alignments[0].hsps:
            top_hit = blast_record.alignments[0]
            hsp = top_hit.hsps[0]
            hit_def = top_hit.hit_def
            identity_pct = (hsp.identities / hsp.align_length) * 100
            genus, genus_species = extract_genus_species(hit_def)

            is_fungal = any(fk.lower() in hit_def.lower() for fk in fungal_keywords)
            included = False

            print(f"  🔎 Top Hit: {hit_def}")
            print(f"  🔗 Identity: {identity_pct:.1f}%")
            print(f"  🍄 Fungal: {'Yes' if is_fungal else 'No'}")

            if is_fungal and identity_pct >= identity_threshold:
                header_name = f"{genus_species} ~{int(identity_pct)}%"
                fasta_out.write(f">{header_name}\n{seq}\n")
                species_counts[genus_species] = species_counts.get(genus_species, 0) + 1
                print(f"  ✅ Included as: {header_name}")
                log_file.write(f"{record.id}\tIncluded\tFungal\t{hit_def} ~{int(identity_pct)}%\n")
            else:
                reason = "Low identity" if is_fungal else "Non-fungal"
                print(f"  ❌ Excluded ({reason})")
                log_file.write(f"{record.id}\tExcluded\t{reason}\t{hit_def} ~{int(identity_pct)}%\n")
        else:
            reason = "No alignment"
            print(f"  ❌ {reason}")
            log_file.write(f"{record.id}\tExcluded\t{reason}\tN/A\n")

        time.sleep(3)

# === Summary ===
print("\n=== ✅ Fungal Matches ===")
for sp, count in species_counts.items():
    print(f"  {sp}: {count}")

print(f"\n📝 Full BLAST log written to: {log_output_path}")
print("🎉 Done.")

