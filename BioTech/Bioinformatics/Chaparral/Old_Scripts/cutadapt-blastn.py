import os
import subprocess
import time
from Bio.Blast import NCBIXML
from Bio import SeqIO
from datetime import datetime

# Platform-specific imports for non-blocking input (optional, not used here)
if os.name == 'nt':
    import msvcrt
else:
    import select
import sys

# === Config ===
input_dir = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
             r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
             r"\fastq_pass\barcode05")
primertype = "FITS"
fwd_primer = "GGAAGTAAAAGTCGTAACAAGG"
rev_primer = "CAAGAGATCCGTTGTTGAAAGTT"
cutadapt_path = r"C:\Users\ketgl\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0\LocalCache\local-packages\Python312\Scripts\cutadapt.exe"

identity_threshold = 75

folder_name = os.path.basename(input_dir.rstrip(os.sep))
input_files = sorted(f for f in os.listdir(input_dir) if f.endswith("_trimmed.fastq"))

if not input_files:
    print("No _trimmed.fastq files found.")
    exit()

final_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_final_trimmed_raw.fastq")

# === Step 1: Run Cutadapt only if final trimmed raw fastq does not exist ===
if os.path.exists(final_output_path):
    print(f"⚡ Found existing trimmed FASTQ: {final_output_path}, skipping Cutadapt step.")
else:
    print("🧹 No existing trimmed FASTQ found, running Cutadapt trimming step...")
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

# === Step 2: Convert trimmed FASTQ to FASTA ===
fasta_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_final_trimmed.fasta")
if os.path.exists(fasta_output_path):
    print(f"⚡ Found existing FASTA: {fasta_output_path}, skipping conversion step.")
else:
    print("🔄 Converting FASTQ to FASTA...")
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

# === Helper functions ===

def extract_genus_species(description):
    parts = description.split()
    if len(parts) >= 2:
        return parts[0], f"{parts[0]} {parts[1]}"
    elif parts:
        return parts[0], parts[0]
    return "unknown", "unknown"

# Local BLAST config
blastn_path = r"C:\Program Files\NCBI\blast-2.16.0+\bin\blastn.exe"
blast_db = r"C:\blast\db\nt"

def blast_with_timeout_local(seq, timeout=300):
    """
    Run local blastn on a sequence with timeout.
    Returns parsed BLAST record or None on timeout/error.
    """
    import tempfile

    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as tmp_fasta:
        tmp_fasta.write(f">query\n{seq}\n")
        query_file = tmp_fasta.name

    with tempfile.NamedTemporaryFile(mode='r', delete=False, suffix=".xml") as tmp_xml:
        output_file = tmp_xml.name

    cmd = [
        blastn_path,
        "-query", query_file,
        "-db", blast_db,
        "-outfmt", "5",  # XML output
        "-max_target_seqs", "1",
        "-out", output_file
    ]

    try:
        subprocess.run(cmd, timeout=timeout, check=True)
    except subprocess.TimeoutExpired:
        print(f"  ⏰ BLAST timed out after {timeout}s. Skipping.")
        cleanup_files(query_file, output_file)
        return None
    except subprocess.CalledProcessError as e:
        print(f"  ❌ BLAST error: {e}")
        cleanup_files(query_file, output_file)
        return None

    try:
        with open(output_file) as result_handle:
            blast_record = NCBIXML.read(result_handle)
    except Exception as e:
        print(f"  ❌ BLAST parse error: {e}")
        blast_record = None

    cleanup_files(query_file, output_file)
    return blast_record

def cleanup_files(*files):
    import os
    for f in files:
        try:
            os.remove(f)
        except Exception:
            pass

# === Filtering criteria for fungi ===
fungal_keywords = [
    "Fungi", "Ascomycota", "Basidiomycota", "Zygomycota", "Chytridiomycota",
    "Glomeromycota", "Microsporidia", "fungus", "mold", "yeast", "mycota", "mycetes"
]

# === BLAST & Output paths ===
fasta_blast_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_fungi_only.fasta")
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_blast_log_{timestamp}.txt")

species_counts = {}

with open(fasta_blast_output_path, "w") as fasta_out, \
     open(log_output_path, "w") as log_file, \
     open(fasta_output_path, "r") as fasta_in:

    log_file.write("ID\tStatus\tReason\tTop_Hit\n")

    for record in SeqIO.parse(fasta_in, "fasta"):
        seq = record.seq
        seq_len = len(seq)
        calculated_timeout = min(480, max(300, int(seq_len / 200)))

        print(f"\n🔬 BLASTing: {record.id} (seq length: {seq_len}, timeout: {calculated_timeout}s)")

        blast_record = blast_with_timeout_local(seq, timeout=calculated_timeout)
        if blast_record is None:
            reason = "Timeout or error during local BLAST"
            print(f"  ⏩ {reason}")
            log_file.write(f"{record.id}\tExcluded\t{reason}\tN/A\n")
            continue

        if blast_record.alignments and blast_record.alignments[0].hsps:
            top_hit = blast_record.alignments[0]
            hsp = top_hit.hsps[0]
            hit_def = top_hit.hit_def
            identity_pct = (hsp.identities / hsp.align_length) * 100
            genus, genus_species = extract_genus_species(hit_def)
            is_fungal = any(fk.lower() in hit_def.lower() for fk in fungal_keywords)

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

        time.sleep(3)  # Respect any policies or avoid overloading

# === Summary ===
print("\n=== ✅ Fungal Matches ===")
for sp, count in species_counts.items():
    print(f"  {sp}: {count}")

print(f"\n📝 Full BLAST log written to: {log_output_path}")
print("🎉 Done.")
