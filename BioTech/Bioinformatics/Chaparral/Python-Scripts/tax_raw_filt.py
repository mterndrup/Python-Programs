import os
import subprocess
from datetime import datetime

# Paths - adjust as needed
base_dir = r"/BioTech/Bioinformatics/Chaparral/Round2"
logs_dir = r"/BioTech/Bioinformatics/Chaparral/Logs"

vsearch_path = r"/BioTech/Bioinformatics/Tools/VSEARCH/vsearch-2.30.0-win-x86_64/bin/vsearch.exe"
sintax_ref_fasta = r"C:/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/Bioinformatics/Tools/UNITE/QIIME/sh_qiime_release_19.02.2025/sh_refs_qiime_ver10_dynamic_19.02.2025_sintax.fasta"


def run_vsearch_sintax(input_fasta, output_path, sintax_db, cutoff=0.95):
    """Run vsearch SINTAX classifier with stderr streaming."""
    cmd = [
        vsearch_path,
        "--sintax", input_fasta,
        "--db", sintax_db,
        "--tabbedout", output_path,
        "--sintax_cutoff", str(cutoff)
    ]

    # Use Popen to capture stderr line-by-line (vsearch prints progress to stderr)
    process = subprocess.Popen(cmd, stderr=subprocess.PIPE, text=True)

    # Print stderr output as it arrives, to avoid buffered/repeated output
    for line in process.stderr:
        print(line, end='')  # vsearch progress lines include newlines

    process.wait()
    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, cmd)


def parse_vsearch_sintax_output(sintax_output_path):
    """Parse vsearch sintax tabbed output to map seq_id -> kingdom."""
    kingdom_map = {}
    with open(sintax_output_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 2:
                # Skip malformed or incomplete lines
                continue
            seq_id = parts[0].split()[0]  # sequence ID (clean)
            taxonomy_str = parts[1]

            kingdom = "Unknown"
            # taxonomy_str example: "d:Fungi(1),p:Ascomycota(0.9),..."
            for entry in taxonomy_str.split(','):
                if entry.startswith('d:'):
                    start = entry.find('(')
                    if start != -1:
                        kingdom = entry[2:start]
                    else:
                        kingdom = entry[2:]
                    break
            kingdom_map[seq_id] = kingdom
    return kingdom_map


def rewrite_fasta_with_kingdom(input_fasta, output_fasta, kingdom_map):
    with open(input_fasta) as fin, open(output_fasta, 'w') as fout:
        header = None
        seq_lines = []

        def write_seq():
            if header is None:
                return
            seq = ''.join(seq_lines)
            clean_id = header.split()[0].lstrip(">")
            kingdom = kingdom_map.get(clean_id, "Unknown")
            new_header = f">{clean_id} kingdom='{kingdom}'"
            fout.write(new_header + "\n")
            fout.write(seq + "\n")

        for line in fin:
            if line.startswith(">"):
                if header is not None:
                    write_seq()
                header = line.strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        write_seq()


def convert_barcode_to_prefix(barcode_folder):
    number = barcode_folder.replace("barcode", "")
    return f"b{number}"


def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    for barcode_folder in sorted(os.listdir(base_dir)):
        barcode_path = os.path.join(base_dir, barcode_folder)
        if not (os.path.isdir(barcode_path) and barcode_folder.startswith("barcode")):
            continue

        sample_prefix = convert_barcode_to_prefix(barcode_folder)
        input_fasta = os.path.join(barcode_path, f"{sample_prefix}_combined_trimmed_R_filt.fasta")
        output_fasta = os.path.join(barcode_path, f"{sample_prefix}_combined_trimmed_R_filt_kingdom.fasta")
        sintax_out = os.path.join(barcode_path, f"{sample_prefix}_vsearch_sintax.txt")

        if not os.path.isfile(input_fasta):
            print(f"File not found: {input_fasta}\n")
            continue

        print(f"Processing: {input_fasta}")

        try:
            run_vsearch_sintax(input_fasta, sintax_out, sintax_ref_fasta, cutoff=0.8)
            kingdom_map = parse_vsearch_sintax_output(sintax_out)
            rewrite_fasta_with_kingdom(input_fasta, output_fasta, kingdom_map)

            counts = {}
            for k in kingdom_map.values():
                counts[k] = counts.get(k, 0) + 1

            print("Kingdom counts:")
            for k, c in sorted(counts.items()):
                print(f"  {k}: {c}")
            print()  # Blank line after each sample

        except subprocess.CalledProcessError as e:
            print(f"vsearch error on {input_fasta}: {e}\n")


if __name__ == "__main__":
    main()
