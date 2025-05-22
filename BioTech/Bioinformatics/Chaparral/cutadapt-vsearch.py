import os
import subprocess
from Bio import SeqIO

# === Hardcoded folder path ===
input_dir = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
             r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
             r"\fastq_pass\barcode07")

# === Primers ===
primertype = "FITS"
fwd_primer = "GGAAGTAAAAGTCGTAACAAGG"
rev_primer = "CAAGAGATCCGTTGTTGAAAGTT"

# === Paths to cutadapt and vsearch ===
cutadapt_path = (r"C:\Users\ketgl\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0"
                 r"\LocalCache\local-packages\Python312\Scripts\cutadapt.exe")

vsearch_path = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                r"\Bioinformatics\Chaparral\vsearch-2.30.0-win-x86_64"
                r"\vsearch-2.30.0-win-x86_64\bin\vsearch.exe")

vsearch_db = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
              r"\Bioinformatics\Chaparral\reference_dbs"
              r"\sh_refs_qiime_ver10_dynamic_19.02.2025.fasta")

# === Taxonomy text file path ===
taxonomy_file = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                 r"\Bioinformatics\Chaparral\reference_dbs"
                 r"\sh_taxonomy_qiime_ver10_dynamic_19.02.2025.txt")

# === File prep ===
folder_name = os.path.basename(input_dir.rstrip(os.sep))
input_files = sorted(f for f in os.listdir(input_dir) if f.endswith("_trimmed.fastq"))
if not input_files:
    print("No _trimmed.fastq files found in the folder.")
    exit()

final_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_final_trimmed_vsearch_raw.fastq")

# === Run cutadapt ===
with open(final_output_path, "w") as final_out:
    for filename in input_files:
        input_path = os.path.join(input_dir, filename)
        temp_output_path = os.path.join(input_dir, f"temp_{filename}_cutadapt.fastq")
        cmd = [
            cutadapt_path, "-g", fwd_primer, "-a", rev_primer,
            "-e", "0.1", "--minimum-length", "60",
            "--overlap", "5",
            "--discard-untrimmed",
            "-o", temp_output_path, input_path
        ]

        print(f"🔪 Running Cutadapt on {filename}...")
        try:
            subprocess.run(cmd, check=True)
            with open(temp_output_path, "r") as temp_in:
                final_out.write(temp_in.read())
            os.remove(temp_output_path)
            print(f"✅ Processed {filename}")
        except subprocess.CalledProcessError as e:
            print(f"❌ Cutadapt error on {filename}: {e}")

print(f"🎉 Final trimmed FASTQ saved to: {final_output_path}")

# === FASTQ → FASTA ===
fasta_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_final_trimmed_vsearch.fasta")
with open(final_output_path, "r") as fastq_in, open(fasta_output_path, "w") as fasta_out:
    for i, line in enumerate(fastq_in):
        if i % 4 == 0:
            fasta_out.write(f">{line[1:]}")
        elif i % 4 == 1:
            fasta_out.write(line)

print(f"📄 FASTA file saved to: {fasta_output_path}")

# === Build SH ID → taxonomy mapping from taxonomy text file ===
print("📖 Loading taxonomy mapping from text file...")
id_to_taxonomy = {}
with open(taxonomy_file, "r") as tf:
    for line in tf:
        if line.strip():
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                sh_id = parts[0]
                tax_string = parts[1]
                tax_parts = tax_string.split(';')
                genus_species = " ".join(tax_parts[-2:]) if len(tax_parts) >= 2 else tax_parts[-1]
                id_to_taxonomy[sh_id] = genus_species

# === Run vsearch taxonomy assignment ===
print("🔍 Starting vsearch taxonomy matching...")

# Cutoffs as you requested
identity_cutoff_confident = 85  # lowered from 90
coverage_cutoff = 80  # 80% coverage cutoff

fasta_vsearch_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_final_trimmed_vsearch_annotated.fasta")
summary_log_path = os.path.join(input_dir, f"{folder_name}_{primertype}_taxonomy_summary.log")

species_counts = {}
unmatched_counts = {}

with open(fasta_vsearch_output_path, "w") as fasta_out:
    for record in SeqIO.parse(fasta_output_path, "fasta"):
        seq = str(record.seq)
        query_id = record.id
        print(f"🔎 Searching: {query_id}")

        temp_query_path = os.path.join(input_dir, "temp_query.fasta")
        with open(temp_query_path, "w") as tq:
            tq.write(f">{query_id}\n{seq}\n")

        vsearch_cmd = [
            vsearch_path,
            "--usearch_global", temp_query_path,
            "--db", vsearch_db,
            "--id", "0.7",
            "--top_hits_only",
            "--blast6out", "temp_vsearch_output.txt",
            "--maxaccepts", "1",
            "--maxhits", "1",
            "--strand", "both"
        ]

        try:
            subprocess.run(vsearch_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"❌ vsearch failed: {e}")
            fasta_out.write(f">uncategorized\n{seq}\n")
            unmatched_counts["uncategorized"] = unmatched_counts.get("uncategorized", 0) + 1
            os.remove(temp_query_path)
            continue

        hit_found = False
        if os.path.exists("temp_vsearch_output.txt"):
            with open("temp_vsearch_output.txt", "r") as vout:
                line = vout.readline().strip()
                if line:
                    cols = line.split("\t")
                    subject_id = cols[1]
                    identity_pct = float(cols[2])
                    align_len = int(cols[3])
                    query_len = len(seq)
                    coverage = (align_len / query_len) * 100 if query_len > 0 else 0

                    taxname = id_to_taxonomy.get(subject_id, None)
                    if taxname and identity_pct >= identity_cutoff_confident and coverage >= coverage_cutoff:
                        header_name = f"{taxname} ~{int(identity_pct)}%"
                        species_counts[taxname] = species_counts.get(taxname, 0) + 1
                    elif identity_pct >= 70:  # borderline matches
                        header_name = f"uncategorized ~{int(identity_pct)}%"
                        unmatched_counts[header_name] = unmatched_counts.get(header_name, 0) + 1
                    else:
                        header_name = "uncategorized"
                        unmatched_counts[header_name] = unmatched_counts.get(header_name, 0) + 1

                    hit_found = True

            os.remove("temp_vsearch_output.txt")

        if not hit_found:
            header_name = "uncategorized"
            unmatched_counts[header_name] = unmatched_counts.get(header_name, 0) + 1

        fasta_out.write(f">{header_name}\n{seq}\n")
        os.remove(temp_query_path)

print(f"📁 vsearch-annotated FASTA saved to: {fasta_vsearch_output_path}")

# === Summary ===
summary_lines = ["\n=== Taxonomy Assignment Summary ===\nMatches:"]
for species, count in species_counts.items():
    summary_lines.append(f"  {species}: {count}")
summary_lines.append("\nUnmatched:")
for label, count in unmatched_counts.items():
    summary_lines.append(f"  {label}: {count}")

summary_text = "\n".join(summary_lines)
print(summary_text)
with open(summary_log_path, "w") as log_file:
    log_file.write(summary_text)

print(f"📝 Summary log saved to: {summary_log_path}")
