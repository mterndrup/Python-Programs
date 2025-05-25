import os
import subprocess
import time
import urllib.error
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Entrez
import threading
import queue
import sys
from datetime import datetime

# Platform-specific imports for non-blocking input
if os.name == 'nt':
    import msvcrt
else:
    import select

# === Prompt for primer type ===
primer_options = {
    "FITS": ("GGAAGTAAAAGTCGTAACAAGG", "CAAGAGATCCGTTGTTGAAAGTT"),
    "PITS": ("ATGCGATACTTGGTGTGAAT", "GACGCTTCTCCAGACTACAAT"),
    "18S": ("GTACACACCGCCCGTC", "TGATCCTTCTGCAGGTTCACCTAC"),
    "16S": ("GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"),
}

print("Select primer type:")
for i, key in enumerate(primer_options.keys(), 1):
    print(f"  {i}. {key}")

while True:
    choice = input("Enter number or primer name (FITS, PITS, 18S, 16S): ").strip().upper()
    if choice in primer_options:
        primertype = choice
        break
    elif choice.isdigit() and 1 <= int(choice) <= len(primer_options):
        primertype = list(primer_options.keys())[int(choice)-1]
        break
    else:
        print("Invalid selection. Please try again.")

fwd_primer, rev_primer = primer_options[primertype]

# === Config ===
input_dir = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
             r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
             r"\fastq_pass\barcode10")

cutadapt_path = r"C:\Users\ketgl\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0\LocalCache\local-packages\Python312\Scripts\cutadapt.exe"

identity_threshold = 75
Entrez.api_key = "b1ee6a94722846e3fe858612221221948607"
Entrez.email = "terndrupm@gmail.com"

folder_name = os.path.basename(input_dir.rstrip(os.sep))
input_files = sorted(f for f in os.listdir(input_dir) if f.endswith("_trimmed.fastq"))

if not input_files:
    print("No _trimmed.fastq files found.")
    exit()

final_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_final_trimmed_raw.fastq")

# === Step 1: Cutadapt Trimming ===
if os.path.exists(final_output_path):
    print(f"⚡ Found existing trimmed FASTQ: {final_output_path}, skipping Cutadapt.")
else:
    print("🧹 No trimmed FASTQ found. Running Cutadapt...")
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
                print(f"❌ Cutadapt error on {filename}: {e}")
    print(f"🎉 Trimmed FASTQ saved: {final_output_path}")

# === Step 2: Convert to FASTA ===
fasta_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_final_trimmed.fasta")
if os.path.exists(fasta_output_path):
    print(f"⚡ Found existing FASTA: {fasta_output_path}, skipping conversion.")
else:
    print("🔄 Converting to FASTA...")
    with open(final_output_path, "r") as fastq_in, open(fasta_output_path, "w") as fasta_out:
        for i, line in enumerate(fastq_in):
            line = line.strip()
            if i % 4 == 0:
                fasta_out.write(f">{line[1:]}\n")
            elif i % 4 == 1:
                fasta_out.write(line + "\n")
    print(f"📂 FASTA saved: {fasta_output_path}")

# === Helper Functions ===
def check_enter_pressed():
    if os.name == 'nt':
        return msvcrt.kbhit() and msvcrt.getwch() in ['\r', '\n']
    else:
        dr, _, _ = select.select([sys.stdin], [], [], 0)
        return dr and sys.stdin.readline().strip() == ""

def blast_with_timeout(seq, timeout=300):
    q = queue.Queue()
    def target():
        try:
            result = NCBIWWW.qblast(
                "blastn", "nt", seq, hitlist_size=1,
                entrez_query=None  # We'll filter later based on taxonomy
            )
            q.put(result)
        except Exception as e:
            q.put(e)
    thread = threading.Thread(target=target)
    thread.daemon = True
    thread.start()
    print(f"  ⏳ BLAST started ({primertype} filter, {timeout}s). Press ENTER to skip...")
    for _ in range(timeout * 10):
        if not thread.is_alive():
            break
        if check_enter_pressed():
            print("  ⏩ Skipped by user.")
            return None
        time.sleep(0.1)
    if thread.is_alive():
        print("  ⏰ Timed out.")
        return None
    result = q.get()
    if isinstance(result, Exception):
        raise result
    return result

def extract_genus_species(description):
    parts = description.split()
    if len(parts) >= 2:
        return parts[0], f"{parts[0]} {parts[1]}"
    elif parts:
        return parts[0], parts[0]
    return "unknown", "unknown"

def taxonomy_filter(hit_def, accession, primertype):
    # Quick check for uncultured fungus for FITS
    if primertype == "FITS" and "uncultured fungus" in hit_def.lower():
        return True
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if not records:
            return False
        lineage = records[0].get("GBSeq_taxonomy", "").lower()
        if primertype == "FITS":
            return "fungi" in lineage
        elif primertype == "PITS":
            return "viridiplantae" in lineage or "plantae" in lineage
        elif primertype == "18S":
            return "eukaryota" in lineage
        elif primertype == "16S":
            return ("bacteria" in lineage) or ("archaea" in lineage)
        else:
            return False
    except Exception as e:
        print(f"  ⚠️ Taxonomy lookup error for {accession}: {e}")
        return False

# === Step 3: BLAST and Filter ===
fasta_blast_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_filtered.fasta")
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
        timeout = min(480, max(300, int(seq_len / 200)))

        print(f"\n🔬 BLASTing {record.id} (len: {seq_len}, timeout: {timeout}s)")
        try:
            result_handle = blast_with_timeout(seq, timeout=timeout)
        except Exception as e:
            reason = f"BLAST error: {e}"
            print(f"  ❌ {reason}")
            log_file.write(f"{record.id}\tExcluded\t{reason}\tN/A\n")
            continue

        if result_handle is None:
            reason = "User skipped or timeout"
            print(f"  ⏩ {reason}")
            log_file.write(f"{record.id}\tExcluded\t{reason}\tN/A\n")
            continue

        try:
            blast_record = NCBIXML.read(result_handle)
        except Exception as e:
            reason = f"Parse error: {e}"
            print(f"  ❌ {reason}")
            log_file.write(f"{record.id}\tExcluded\t{reason}\tN/A\n")
            continue

        if blast_record.alignments and blast_record.alignments[0].hsps:
            top_hit = blast_record.alignments[0]
            hsp = top_hit.hsps[0]
            hit_def = top_hit.hit_def
            accession = top_hit.accession
            identity_pct = (hsp.identities / hsp.align_length) * 100
            genus, genus_species = extract_genus_species(hit_def)

            print(f"  🔎 Top Hit: {hit_def}")
            print(f"  🔗 Identity: {identity_pct:.1f}%")

            passes_filter = taxonomy_filter(hit_def, accession, primertype)

            if passes_filter and identity_pct >= identity_threshold:
                header = f"{genus_species} ~{int(identity_pct)}%"
                fasta_out.write(f">{header}\n{seq}\n")
                species_counts[genus_species] = species_counts.get(genus_species, 0) + 1
                print(f"  ✅ Included as: {header}")
                log_file.write(f"{record.id}\tIncluded\t{primertype} filter\t{hit_def} ~{int(identity_pct)}%\n")
            else:
                reason = "Low identity" if passes_filter else "Taxonomy filter failed"
                print(f"  ❌ Excluded ({reason})")
                log_file.write(f"{record.id}\tExcluded\t{reason}\t{hit_def} ~{int(identity_pct)}%\n")
        else:
            reason = "No alignment"
            print(f"  ❌ {reason}")
            log_file.write(f"{record.id}\tExcluded\t{reason}\tN/A\n")

        time.sleep(3)  # NCBI rate limit

# === Summary ===
print("\n=== ✅ Filtered Matches ===")
for sp, count in species_counts.items():
    print(f"  {sp}: {count}")

print(f"\n📝 BLAST log saved to: {log_output_path}")
print("🎉 Done.")
