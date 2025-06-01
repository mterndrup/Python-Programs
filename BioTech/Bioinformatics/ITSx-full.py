import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from datetime import datetime

def get_base_dir(project_choice):
    if project_choice == "1":
        return "/mnt/c/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/Bioinformatics/Chaparral/Wildfire/fastq_pass"
    elif project_choice == "2":
        return "/mnt/c/Users/ketgl/OneDrive/Documents/Glowing-Fungi"
    elif project_choice == "3":
        return "/mnt/c/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/Bioinformatics/Plants/Cacti"
    elif project_choice == "4":  # Added Aya
        return "/mnt/c/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/Bioinformatics/Plants/Aya/DNA_Sequences"
    else:
        return None

def get_fasta_path(project_choice, base_dir, barcode_num, mode, use_18s=False):
    if project_choice == "1" or project_choice == "4":
        # Chaparral or Aya: barcode folder and fasta inside it
        barcode_folder = f"barcode{barcode_num}"
        folder_path = os.path.join(base_dir, barcode_folder)
        if mode == "plant" and use_18s:
            # Look for _PITS_18S_combo.fasta file in barcode folder
            for file in os.listdir(folder_path):
                if file.endswith("_PITS_18S_combo.fasta"):
                    return os.path.join(folder_path, file)
            print("Warning: _PITS_18S_combo.fasta file not found, using default naming")
            return os.path.join(folder_path, f"barcode{barcode_num}_PITS_18S_combo.fasta")
        else:
            suffix = "PITS" if mode == "plant" else "FITS"
            fasta_file = f"barcode{barcode_num}_{suffix}.fasta"
            return os.path.join(folder_path, fasta_file)
    elif project_choice in ("2", "3"):
        return os.path.join(base_dir, "assembly.fasta")
    else:
        return None

def make_headers_unique(input_fasta, output_fasta):
    print("Making sequence headers unique...")
    seen_ids = set()
    unique_records = []
    for i, record in enumerate(SeqIO.parse(input_fasta, "fasta"), 1):
        new_id = record.id
        if new_id in seen_ids:
            new_id = f"{record.id}_{i}"
        seen_ids.add(new_id)
        record.id = new_id
        record.name = new_id
        record.description = ""
        unique_records.append(record)
    SeqIO.write(unique_records, output_fasta, "fasta")
    print(f"Unique-header FASTA written to: {output_fasta}")
    return output_fasta

def split_long_contigs(input_fasta, max_len=99000):
    base_dir = os.path.dirname(input_fasta)
    base_name = os.path.basename(input_fasta)
    split_fasta = os.path.join(base_dir, base_name.replace(".fasta", "_split.fasta"))

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

def run_itsx(input_fasta, mode, output_dir, is_full_combo=False):
    itsx_path = "/mnt/c/Users/ketgl/Downloads/ITSx_1.1.3/ITSx_1.1.3/ITSx"
    cpus = "4"

    kingdom = "F" if mode == "fungi" else "T"

    # Build output prefix inside output_dir with unique timestamp
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    base_name = os.path.basename(input_fasta).replace(".fasta", "")
    output_prefix = os.path.join(output_dir, f"{base_name}_ITSx_{timestamp}")

    cmd = [
        "perl",
        itsx_path,
        "-i", input_fasta,
        "-o", output_prefix,
        "-t", kingdom,
        "--cpu", cpus,
        "--preserve", "T",
        "--allow_partial",
        "--save_regions", "all"
    ]

    print(f"Running ITSx on {input_fasta} (kingdom={kingdom})")
    try:
        subprocess.run(cmd, check=True)
        print(f"ITSx completed successfully. Output files prefix: {output_prefix}")
    except subprocess.CalledProcessError as e:
        print(f"Error running ITSx: {e}")
        return False

    return output_prefix

def check_full_files_and_rerun(itsx_output_dir, mode, barcode_folder, barcode_num, use_18s, rerun_allowed):
    print(f"Checking .full.fasta files in ITSx output folder: {itsx_output_dir}")
    full_files = [f for f in os.listdir(itsx_output_dir) if f.endswith(".full.fasta")]

    need_rerun = False
    for full_file in full_files:
        full_path = os.path.join(itsx_output_dir, full_file)
        size = os.path.getsize(full_path)
        if size == 0:
            print(f"Found empty .full.fasta file: {full_file}")
            need_rerun = True
            break
        else:
            # Check if fasta file actually contains sequences
            if not any(SeqIO.parse(full_path, "fasta")):
                print(f".full.fasta file has no sequences: {full_file}")
                need_rerun = True
                break

    if need_rerun and rerun_allowed:
        print("Rerunning ITSx on full combo fasta due to empty .full.fasta file...")

        full_combo_fasta = os.path.join(barcode_folder, f"barcode{barcode_num}_full_combo.fasta")
        if not os.path.isfile(full_combo_fasta):
            print(f"Full combo fasta not found: {full_combo_fasta}. Cannot rerun ITSx.")
            return

        unique_full_combo = os.path.join(barcode_folder, f"barcode{barcode_num}_full_combo_combined.fasta")
        make_headers_unique(full_combo_fasta, unique_full_combo)
        split_full_combo = split_long_contigs(unique_full_combo)

        full_combo_output_dir = os.path.join(itsx_output_dir, f"ITSx_fullcombo_{datetime.now().strftime('%Y%m%d%H%M%S')}")
        os.makedirs(full_combo_output_dir, exist_ok=True)

        run_itsx(split_full_combo, mode, full_combo_output_dir, is_full_combo=True)
    elif need_rerun and not rerun_allowed:
        print("ITSx has already been rerun on full combo fasta once. Stopping to avoid infinite loop.")
    else:
        print("All .full.fasta files contain sequences. No rerun needed.")

def main():
    while True:
        print("Select project:")
        print("1) Chaparral")
        print("2) Glowing Mushroom")
        print("3) Cacti")
        print("4) Aya")  # Added Aya option
        project_choice = input("Enter project number (1/2/3/4): ").strip()
        if project_choice in ("1", "2", "3", "4"):
            break
        print("Invalid input. Please enter 1, 2, 3 or 4.")

    barcode_num = None
    if project_choice in ("1", "4"):  # Chaparral and Aya need barcode input
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

    use_18s = False
    if mode == "plant":
        while True:
            response = input("Include 18S sequences? (yes/no): ").strip().lower()
            if response in ("yes", "no"):
                use_18s = (response == "yes")
                break
            print("Please answer 'yes' or 'no'.")

    base_dir = get_base_dir(project_choice)
    if base_dir is None:
        print("Invalid project choice or base directory not set. Exiting.")
        return

    fasta_path = get_fasta_path(project_choice, base_dir, barcode_num if barcode_num else "", mode, use_18s)
    print(f"\nUsing FASTA path:\n{fasta_path}")

    if not fasta_path or not os.path.isfile(fasta_path):
        print(f"FASTA file not found:\n{fasta_path}\nExiting.")
        return

    barcode_folder = os.path.dirname(fasta_path)

    # Make unique headers fasta
    base_name = os.path.basename(fasta_path).replace(".fasta", "")
    unique_headers_fasta = os.path.join(barcode_folder, f"{base_name}_combined.fasta")
    make_headers_unique(fasta_path, unique_headers_fasta)

    # Split long contigs fasta
    split_fasta = split_long_contigs(unique_headers_fasta)

    # Create ITSx output folder inside base_dir
    itsx_dir = os.path.join(barcode_folder, "ITSx")
    os.makedirs(itsx_dir, exist_ok=True)

    # Run ITSx first time
    first_output_subdir = os.path.join(itsx_dir, f"ITSx_run_{datetime.now().strftime('%Y%m%d%H%M%S')}")
    os.makedirs(first_output_subdir, exist_ok=True)

    output_prefix = run_itsx(split_fasta, mode, first_output_subdir, is_full_combo=False)
    if not output_prefix:
        print("ITSx failed, exiting.")
        return

    # Check .full.fasta files and rerun if empty and rerun allowed (only once)
    check_full_files_and_rerun(first_output_subdir, mode, barcode_folder, barcode_num if barcode_num else "", use_18s, rerun_allowed=True)

if __name__ == "__main__":
    main()
