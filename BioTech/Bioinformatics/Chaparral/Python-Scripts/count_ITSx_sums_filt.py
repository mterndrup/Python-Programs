import os
import csv
from datetime import datetime

BASE_FOLDER = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\Round2"
LOG_FOLDER = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\Logs"

ITS_REGIONS = ['5_8S', 'ITS1', 'ITS2', 'LSU', 'SSU', 'full']

def read_fasta_to_dict(fasta_path):
    """Read fasta sequences into dict: {seqID: sequence}"""
    seq_dict = {}
    with open(fasta_path, 'r', encoding='utf-8') as f:
        seq_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id:
                    seq_dict[seq_id] = ''.join(seq_lines)
                # seqID is header up to first space
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id:
            seq_dict[seq_id] = ''.join(seq_lines)
    return seq_dict

def count_csv_rows(csv_path):
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        return sum(1 for _ in reader)

def count_its_regions(csv_path):
    """Count sum of 1's for each ITS region in filtered csv."""
    counts = {region: 0 for region in ITS_REGIONS}
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for region in ITS_REGIONS:
                counts[region] += int(row.get(region, '0'))
    return counts

def extract_its_regions(barcode_folder, barcode_name):
    """
    Extract sequences for ITS regions indicated as 1 in filtered csv
    from ITSx fasta files in ITSx-v1 folder.
    Save results to CSV with seqID and ITS sequences columns.
    """
    print(f"  Extracting ITS sequences for {barcode_name}...")
    filtered_csv = os.path.join(barcode_folder, f"{barcode_name}_ITSx_seq_region_filtered.csv")
    itsx_folder = os.path.join(barcode_folder, "ITSx-v1")

    # Load filtered CSV data
    filtered_rows = []
    with open(filtered_csv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            filtered_rows.append(row)

    # Load ITSx fasta files into dicts per region
    fasta_files = {}
    for region in ITS_REGIONS:
        fasta_name = f"{barcode_name.replace('barcode', 'b')}_combined_trimmed_R_filt_ITSx.{region}.fasta"
        fasta_path = os.path.join(itsx_folder, fasta_name)
        if os.path.isfile(fasta_path):
            print(f"    Loading {region} from {fasta_path}")
            fasta_files[region] = read_fasta_to_dict(fasta_path)
        else:
            print(f"    WARNING: Missing ITSx fasta file for {region} at {fasta_path}")

    # Prepare output CSV rows
    output_rows = []
    header = ['seqID'] + ITS_REGIONS
    for row in filtered_rows:
        seq_id = row['seqID']
        out_row = {'seqID': seq_id}
        for region in ITS_REGIONS:
            if row.get(region, '0') == '1':
                seq = fasta_files.get(region, {}).get(seq_id, '')
                out_row[region] = seq
            else:
                out_row[region] = ''
        output_rows.append(out_row)

    # Write output CSV
    out_csv = os.path.join(barcode_folder, f"{barcode_name}_ITS_regions_extracted.csv")
    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"    Extracted ITS regions saved to: {out_csv}")

def process_barcode_folder(barcode_folder):
    barcode_name = os.path.basename(barcode_folder)
    print(f"[PROCESSING] {barcode_name}")

    # Define expected files
    # original fasta - not used for filtering summary now, so no counts done on this
    original_fasta = os.path.join(barcode_folder, f"{barcode_name.replace('barcode','b')}_combined_trimmed_R_filt.fasta")
    filtered_fasta = os.path.join(barcode_folder, f"{barcode_name.replace('barcode','b')}_combined_trimmed_R_filt_kingdom.fasta")
    original_csv = os.path.join(barcode_folder, f"{barcode_name}_ITSx_seq_region_counts.csv")
    filtered_csv = os.path.join(barcode_folder, f"{barcode_name}_ITSx_seq_region_filtered.csv")

    # Check required files exist
    missing = False
    for file_path in [original_csv, filtered_csv]:
        if not os.path.isfile(file_path):
            print(f"  ERROR: Missing file {file_path}, skipping {barcode_name}")
            missing = True
    if missing:
        return None

    # Count rows in CSVs
    total_csv_orig = count_csv_rows(original_csv)
    total_csv_filt = count_csv_rows(filtered_csv)
    filtered_out_csv = total_csv_orig - total_csv_filt
    filtered_out_pct = (filtered_out_csv / total_csv_orig * 100) if total_csv_orig else 0.0

    # Count ITS regions in filtered csv
    its_region_counts = count_its_regions(filtered_csv)

    # Extract ITS sequences
    extract_its_regions(barcode_folder, barcode_name)

    # Write summary file
    summary_path = os.path.join(barcode_folder, f"{barcode_name}_summary.txt")
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write(f"Summary for {barcode_name}\n\n")
        f.write("CSV Filtering:\n")
        f.write(f"Original CSV rows: {total_csv_orig}\n")
        f.write(f"Filtered CSV rows: {total_csv_filt}\n")
        f.write(f"Filtered out rows: {filtered_out_csv} ({filtered_out_pct:.2f}%)\n\n")

        f.write("ITS Region Counts (filtered CSV):\n")
        for region in ITS_REGIONS:
            f.write(f"{region}: {its_region_counts.get(region, 0)}\n")

        f.write("\nNote: This filtering was to keep fungi sequences only.\n")

    print(f"  Summary saved to: {summary_path}")

    return {
        'barcode': barcode_name,
        'original_csv_rows': total_csv_orig,
        'filtered_csv_rows': total_csv_filt,
        'filtered_out_rows': filtered_out_csv,
        'filtered_out_pct': filtered_out_pct,
        'its_region_counts': its_region_counts
    }

def write_global_log(all_results):
    os.makedirs(LOG_FOLDER, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(LOG_FOLDER, f"global_summary_{timestamp}.txt")

    # Calculate combined totals
    total_orig = sum(r['original_csv_rows'] for r in all_results)
    total_filt = sum(r['filtered_csv_rows'] for r in all_results)
    total_filtered_out = sum(r['filtered_out_rows'] for r in all_results)
    filtered_out_pct = (total_filtered_out / total_orig * 100) if total_orig else 0.0

    # Sum ITS region counts
    combined_its_counts = {region: 0 for region in ITS_REGIONS}
    for r in all_results:
        for region in ITS_REGIONS:
            combined_its_counts[region] += r['its_region_counts'].get(region, 0)

    lines = [
        f"Global Summary for all barcodes\n",
        f"Total Original CSV rows: {total_orig}",
        f"Total Filtered CSV rows: {total_filt}",
        f"Total Filtered out rows: {total_filtered_out} ({filtered_out_pct:.2f}%)\n",
        "Combined ITS Region Counts (filtered CSV):"
    ]
    for region in ITS_REGIONS:
        lines.append(f"{region}: {combined_its_counts[region]}")
    lines.append("\nNote: This filtering was to keep fungi sequences only.\n")

    with open(log_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))

    print(f"[GLOBAL SUMMARY] saved to: {log_path}")

def main():
    all_results = []
    for entry in os.listdir(BASE_FOLDER):
        barcode_path = os.path.join(BASE_FOLDER, entry)
        if os.path.isdir(barcode_path) and entry.startswith("barcode"):
            result = process_barcode_folder(barcode_path)
            if result:
                all_results.append(result)

    if all_results:
        write_global_log(all_results)
    else:
        print("No barcodes processed successfully.")

if __name__ == "__main__":
    main()
