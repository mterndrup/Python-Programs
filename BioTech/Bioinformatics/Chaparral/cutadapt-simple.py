import os
import subprocess
import sys

# --- Config ---
fastq_base_dir = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                  r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
                  r"\fastq_pass")

cutadapt_path = r"C:\Users\ketgl\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0\LocalCache\local-packages\Python312\Scripts\cutadapt.exe"

primers = {
    "FITS": {
        "fwd": "GGAAGTAAAAGTCGTAACAAGG",
        "rev": "CAAGAGATCCGTTGTTGAAAGTT"
    },
    # "PITS": {
    #     "fwd": "ATGCGATACTTGGTGTGAAT",
    #     "rev": "GACGCTTCTCCAGACTACAAT"
    # },
    # "18S": {
    #     "fwd": "GTACACACCGCCCGTC",
    #     "rev": "TGATCCTTCTGCAGGTTCACCTAC"
    # }
}

# === Helper Functions ===

def run_cutadapt(input_path, primer_name, fwd_primer, rev_primer, out_dir):
    basename = os.path.splitext(os.path.basename(input_path))[0]
    output_path = os.path.join(out_dir, f"{basename}_{primer_name}_cutadapt.fastq")
    cmd = [
        cutadapt_path,
        "-g", f"{fwd_primer}",
        "-a", f"{rev_primer}",
        "-e", "0.1",
        "--minimum-length", "60",
        "--overlap", "7",
        #"--discard-untrimmed",
        "-o", output_path,
        input_path
    ]
    print(f"🔪 Running Cutadapt for {basename} with {primer_name} primers")
    subprocess.run(cmd, check=True)
    return output_path

def get_suffix(filename):
    name = os.path.splitext(filename)[0]
    parts = name.split('_')
    suffix_parts = parts[-4:]
    suffix = "_" + "_".join(suffix_parts)
    return suffix

def fastq_to_fasta_records(fastq_path, primer_name, input_file, min_length=100):
    records = []
    try:
        suffix = get_suffix(input_file)
        with open(fastq_path, "r") as fin:
            lines = fin.readlines()
            for i in range(0, len(lines), 4):
                header_line = lines[i].strip()
                sequence = lines[i + 1].strip()
                if len(sequence) < min_length:
                    continue  # Skip sequences shorter than min_length
                unique_id = header_line.split()[0][1:]  # remove '@'
                full_id = f"{unique_id}_x{suffix}"
                records.append(f">{full_id}\n{sequence}\n")
    except Exception as e:
        print(f"⚠️ Warning: Failed reading {fastq_path}: {e}")
    return records


# === MAIN PROCESS ===

barcode_dirs = [d for d in os.listdir(fastq_base_dir) if d.startswith("barcode") and os.path.isdir(os.path.join(fastq_base_dir, d))]

if not barcode_dirs:
    print("🚫 No barcode folders found.")
    sys.exit(1)

for barcode_dir in barcode_dirs:
    barcode_path = os.path.join(fastq_base_dir, barcode_dir, "Raw")
    if not os.path.exists(barcode_path):
        print(f"📁 Skipping {barcode_dir} (no Raw folder)")
        continue

    input_files = sorted(f for f in os.listdir(barcode_path) if f.endswith("_trimmed_filtered.fastq"))
    if not input_files:
        print(f"📭 No _trimmed.fastq files found in {barcode_dir}")
        continue

    all_fasta_records = []
    print(f"\n🔍 Processing {barcode_dir}")

    for input_file in input_files:
        input_path = os.path.join(barcode_path, input_file)

        for primer_name, primers_seq in primers.items():
            cutadapt_out = run_cutadapt(input_path, primer_name, primers_seq["fwd"], primers_seq["rev"], barcode_path)
            records = fastq_to_fasta_records(cutadapt_out, primer_name, input_file, min_length=100)
            print(f"📥 {len(records)} sequences added from {primer_name} ({input_file})")
            all_fasta_records.extend(records)
            os.remove(cutadapt_out)

    # Write per-barcode FASTA inside barcode folder but outside Raw
    barcode_output_path = os.path.join(fastq_base_dir, barcode_dir, f"{barcode_dir}_FITS_full.fasta")
    with open(barcode_output_path, "w") as fout:
        fout.writelines(all_fasta_records)

    print(f"✅ Final FASTA for {barcode_dir} written: {barcode_output_path}")
