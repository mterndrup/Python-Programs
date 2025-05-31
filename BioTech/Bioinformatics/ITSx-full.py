import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_base_dir(project_choice):
    if project_choice == "1":
        return "/mnt/c/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/Bioinformatics/Chaparral/wildfirePlants-DNA-nanopore-sequence/fastq_pass"
    elif project_choice == "2":
        return "/mnt/c/Users/ketgl/OneDrive/Documents/Glowing-Fungi"
    elif project_choice == "3":
        return "/mnt/c/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/Bioinformatics/Plants/Cacti"
    else:
        return None

def get_fasta_path(project_choice, base_dir, barcode_num, mode):
    if project_choice == "1":
        suffix = "PITS" if mode == "plant" else "FITS"
        barcode_folder = f"barcode{barcode_num}"
        fasta_file = f"barcode{barcode_num}_{suffix}.fasta"
        return os.path.join(base_dir, barcode_folder, fasta_file)
    elif project_choice in ("2", "3"):
        return os.path.join(base_dir, "assembly.fasta")
    else:
        return None

def split_long_contigs(input_fasta, max_len=99000):
    base_dir = os.path.dirname(input_fasta)
    split_fasta = os.path.join(base_dir, "assembly_contigs_split.fasta")

    print("Splitting long contigs if needed...")
    split_records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = record.seq
        if len(seq) <= max_len:
            split_records.append(record)
        else:
            for i in range(0, len(seq), max_len):
                chunk = seq[i:i + max_len]
                chunk_id = f"{record.id}_part_{i // max_len + 1}"
                split_records.append(SeqRecord(chunk, id=chunk_id, description=""))

    SeqIO.write(split_records, split_fasta, "fasta")
    print(f"Split contigs written to: {split_fasta}")
    return split_fasta

def run_itsx(split_fasta, mode, base_dir):
    itsx_path = "/mnt/c/Users/ketgl/Downloads/ITSx_1.1.3/ITSx_1.1.3/ITSx"
    cpus = "4"

    suffix = "FITS" if mode == "fungi" else "PITS"
    itsx_output_dir = os.path.join(base_dir, "ITSx")
    os.makedirs(itsx_output_dir, exist_ok=True)
    output_prefix = os.path.join(itsx_output_dir, f"ITSx_out_{suffix}")

    kingdom = "Fungi" if mode == "fungi" else "Plantae"

    cmd = [
        "perl",
        itsx_path,
        "-i", split_fasta,
        "-o", output_prefix,
        "--cpu", cpus,
        "--preserve", "T",
        "--allow_partial",
        "--only_full",
        "--kingdom", kingdom
    ]

    print(f"Running ITSx restricted to kingdom: {kingdom}")
    try:
        subprocess.run(cmd, check=True)
        print(f"ITSx completed successfully. Output: {output_prefix}*")
    except subprocess.CalledProcessError as e:
        print(f"Error running ITSx: {e}")

def main():
    while True:
        print("Select project:")
        print("1) Chaparral")
        print("2) Glowing Mushroom")
        print("3) Cacti")
        project_choice = input("Enter project number (1/2/3): ").strip()
        if project_choice in ("1", "2", "3"):
            break
        print("Invalid input. Please enter 1, 2, or 3.")

    barcode_num = None
    if project_choice == "1":
        while True:
            barcode_num = input("Enter barcode number (e.g., 1 or 05): ").strip()
            if barcode_num.isdigit():
                barcode_num = barcode_num.zfill(2)
                break
            print("Invalid barcode number. Please enter digits only.")

    while True:
        mode = input("Run ITSx in plant or fungi mode? (plant/fungi): ").strip().lower()
        if mode in ("plant", "fungi"):
            break
        print("Invalid input. Please type 'plant' or 'fungi'.")

    base_dir = get_base_dir(project_choice)
    fasta_path = get_fasta_path(project_choice, base_dir, barcode_num if barcode_num else "", mode)

    print(f"\nUsing FASTA path:\n{fasta_path}")

    if not fasta_path or not os.path.isfile(fasta_path):
        print(f"FASTA file not found:\n{fasta_path}\nExiting.")
        return

    split_fasta = split_long_contigs(fasta_path)
    run_itsx(split_fasta, mode, os.path.dirname(split_fasta))

if __name__ == "__main__":
    main()
