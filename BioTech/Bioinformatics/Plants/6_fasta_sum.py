import re
from pathlib import Path

# === Config ===
BASE_DIR = Path(r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus-All")

# Strain lookup dictionary
strain_lookup = {
    "Lumberjack": 33, "Sharxx": 36, "Reynolds": 41, "Misplant Crested": 42,
    "Kazumi": 44, "Eileen": 47, "SLO": 34, "Zelly 3": 38, "JVD": 39,
    "Psycho 0": 40, "Apple Crested": 46, "Buena Creek (BC)": 45, "Mimis Pachanoi Var": 48,
    "Floyd Old": 86, "Kate Old": 85, "Oceanside Old": 84, "Oceanside": 73,
    "Kate": 74, "Floyd": 75, "Morrison": 77, "Zelly Boobie": 76, "SS01": 78,
    "Daniel's": 79, "Penis Cactus": 80, "Eve's Pin": 83
}

barcode_to_strain = {v: k for k, v in strain_lookup.items()}


def parse_fasta(file_path):
    records = []
    header, seq_lines = None, []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq_lines)))
                header, seq_lines = line, []
            else:
                seq_lines.append(line)
        if header:
            records.append((header, "".join(seq_lines)))
    return records


def get_seq_id(header):
    match = re.match(r">([0-9a-fA-F]+)", header)
    return match.group(1) if match else "unknown"


def get_top_records(barcode_folder, region, top_n):
    logs_dir = barcode_folder / "Logs"
    if not logs_dir.exists():
        return []
    fasta_files = list(logs_dir.glob(f"*_{region}_clusters_selected.fasta"))
    if not fasta_files:
        return []
    records = parse_fasta(fasta_files[0])
    return records[:top_n]


def main():
    # --- Ask user for regions ---
    print("Select region to extract:")
    print("1. ITS2")
    print("2. LSU")
    print("3. Both")
    choice = input("Enter 1, 2, or 3: ").strip()

    if choice == "1":
        regions = ["ITS2"]
    elif choice == "2":
        regions = ["LSU"]
    elif choice == "3":
        regions = ["ITS2", "LSU"]
    else:
        print("Invalid choice. Exiting.")
        return

    # --- Number of sequences per barcode ---
    try:
        top_n = int(input("How many sequences per barcode? ").strip())
    except ValueError:
        print("Invalid number")
        return

    # --- Concatenate results option ---
    concat_choice = input("Concatenate results into one FASTA? (y/n): ").strip().lower()
    concatenate = concat_choice == "y"

    # --- Barcode selection ---
    print("\nDo you want to process:")
    print("1. All barcodes")
    print("2. A single barcode")
    bc_choice = input("Enter 1 or 2: ").strip()

    barcode_dirs = []
    if bc_choice == "1":
        # all barcodes
        barcode_dirs = [d for d in BASE_DIR.iterdir() if d.is_dir() and d.name.startswith("barcode")]
    elif bc_choice == "2":
        try:
            bc_num = int(input("Enter the barcode number (e.g. 33): ").strip())
            target_dir = BASE_DIR / f"barcode{bc_num}"
            if target_dir.exists():
                barcode_dirs = [target_dir]
            else:
                print(f"Barcode {bc_num} not found in {BASE_DIR}")
                return
        except ValueError:
            print("Invalid barcode number.")
            return
    else:
        print("Invalid choice.")
        return

    combined_output_lines = []

    for barcode_dir in barcode_dirs:
        barcode_num = int("".join(filter(str.isdigit, barcode_dir.name)))
        strain_name = barcode_to_strain.get(barcode_num, f"Unknown{barcode_num}")

        its_records = get_top_records(barcode_dir, "ITS2", top_n) if "ITS2" in regions else []
        lsu_records = get_top_records(barcode_dir, "LSU", top_n) if "LSU" in regions else []

        if "ITS2" in regions and "LSU" in regions and concatenate:
            # Combine ITS + LSU per index
            for i in range(max(len(its_records), len(lsu_records))):
                its_header, its_seq = its_records[i] if i < len(its_records) else ("", "")
                lsu_header, lsu_seq = lsu_records[i] if i < len(lsu_records) else ("", "")

                its_id = get_seq_id(its_header) if its_header else "unknown"
                lsu_id = get_seq_id(lsu_header) if lsu_header else "unknown"

                combined_header = f">{strain_name.replace(' ', '')}({barcode_num})_{its_id}_ITS_{lsu_id}_LSU"
                combined_seq = its_seq + lsu_seq

                combined_output_lines.append(combined_header)
                combined_output_lines.append(combined_seq)
        else:
            # Either single region or not concatenating
            for region_name, records in zip(regions, [its_records, lsu_records]):
                for header, seq in records:
                    seq_id = get_seq_id(header)
                    output_header = f">{strain_name.replace(' ', '')}({barcode_num})_{seq_id}_{region_name if region_name=='LSU' else 'ITS'}"
                    combined_output_lines.append(output_header)
                    combined_output_lines.append(seq)

        print(f"✔ Processed barcode {barcode_dir.name} ({strain_name})")

    # --- Write output ---
    if concatenate or len(regions) == 2:
        output_file = BASE_DIR / f"Cactus_top{top_n}_ITS_LSU_combined.fasta"
    else:
        output_file = BASE_DIR / f"Cactus_top{top_n}_{'_'.join(regions)}_summary.fasta"

    with open(output_file, "w") as out:
        out.write("\n".join(combined_output_lines))

    print(f"\n✅ Wrote FASTA: {output_file}")


if __name__ == "__main__":
    main()
