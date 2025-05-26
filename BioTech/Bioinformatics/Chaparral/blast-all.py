import os
import time
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Entrez
import threading
import sys
from datetime import datetime
import csv
import glob
import re

identity_threshold = 85
Entrez.api_key = "b1ee6a94722846e3fe858612221221948607"
Entrez.email = "terndrupm@gmail.com"
primertype = "FITS"

# Base path where barcode folders reside
base_dir = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence\fastq_pass"

def get_barcode_folders(base_dir):
    # List all folders matching barcodeXX pattern
    folders = []
    pattern = re.compile(r"barcode(\d{2})$")
    for entry in os.listdir(base_dir):
        full_path = os.path.join(base_dir, entry)
        if os.path.isdir(full_path):
            m = pattern.match(entry)
            if m:
                folders.append((int(m.group(1)), full_path))
    folders.sort(key=lambda x: x[0])
    return folders

def run_blast_for_barcode(barcode_num, input_dir):
    barcode = str(barcode_num).zfill(2)

    print(f"\n=== Processing barcode{barcode} ===")

    input_fasta_path = os.path.join(input_dir, f"barcode{barcode}_{primertype}.fasta")
    if not os.path.exists(input_fasta_path):
        print(f"❌ Input FASTA file not found: {input_fasta_path}. Skipping barcode{barcode}.")
        return True  # Skip and continue to next barcode

    folder_name = os.path.basename(input_dir.rstrip(os.sep))

    fasta_blast_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_BLAST.fasta")

    log_files = sorted(glob.glob(os.path.join(input_dir, f"{folder_name}_{primertype}_blast_log_*.txt")))

    if log_files:
        log_output_path = log_files[-1]
        print(f"⚡ Resuming from existing log: {log_output_path}")
    else:
        log_output_path = os.path.join(input_dir, f"{folder_name}_{primertype}_blast_log_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.txt")
        print(f"🆕 Starting new log file: {log_output_path}")

    summary_csv_path = os.path.join(input_dir, f"{folder_name}_sum.csv")

    species_counts = {}
    processed_ids = set()

    if os.path.exists(fasta_blast_output_path):
        for record in SeqIO.parse(fasta_blast_output_path, "fasta"):
            processed_ids.add(record.id)

    if os.path.exists(log_output_path):
        with open(log_output_path, "r", encoding="utf-8") as log_file:
            header_line = log_file.readline()
            for line in log_file:
                if line.strip():
                    rec_id = line.split("\t")[0]
                    processed_ids.add(rec_id)

    if os.path.exists(summary_csv_path):
        with open(summary_csv_path, "r", encoding="utf-8") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                species_counts[row["Species"]] = int(row["Count"])

    def blast_with_timeout(seq, timeout=300):
        blast_result = [None]

        def blast_task():
            try:
                print("🔁 Starting BLAST request...")
                result_handle = NCBIWWW.qblast("blastn", "nt", seq, hitlist_size=1,
                                               entrez_query="txid4751[ORGN]")
                return result_handle
            except Exception as e:
                print(f"❌ Exception during BLAST: {e}")
                return e

        thread = threading.Thread(target=lambda: blast_result.__setitem__(0, blast_task()))
        thread.start()

        start_time = time.time()
        while thread.is_alive():
            elapsed = time.time() - start_time
            if elapsed > timeout:
                print(f"⏰ Timeout after {timeout} seconds (sequence length: {len(seq)}). Marking as skipped.")
                return None
            time.sleep(0.1)

        thread.join(timeout=1)
        result = blast_result[0]

        if result is None:
            print("⚠️ No result returned from BLAST thread — possible crash or premature exit.")
            return None

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

    def is_fungi_kingdom_or_uncultured(hit_def, accession):
        if "uncultured fungus" in hit_def.lower():
            return True
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            if not records:
                return False
            lineage = records[0].get("GBSeq_taxonomy", "")
            return "Fungi" in lineage
        except Exception as e:
            print(f"  ⚠️ Taxonomy lookup error for {accession}: {e}")
            return False

    with open(fasta_blast_output_path, "a") as fasta_out, \
         open(log_output_path, "a", encoding="utf-8") as log_file:

        if os.path.getsize(log_output_path) == 0:
            log_file.write("ID\tStatus\tReason\tTop_Hit\n")

        for record in SeqIO.parse(input_fasta_path, "fasta"):
            if record.id in processed_ids:
                print(f"⏭️ Skipping already-processed: {record.id}")
                continue

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

                is_fungal = is_fungi_kingdom_or_uncultured(hit_def, accession)

                if is_fungal and identity_pct >= identity_threshold:
                    header = f"{genus_species} ~{int(identity_pct)}%"
                    fasta_out.write(f">{record.id}\n{seq}\n")
                    species_counts[genus_species] = species_counts.get(genus_species, 0) + 1
                    print(f"  ✅ Included as: {header}")
                    log_file.write(f"{record.id}\tIncluded\tFungal\t{hit_def} ~{int(identity_pct)}%\n")
                else:
                    reason = "Low identity" if is_fungal else "Non-fungal"
                    print(f"  ❌ Excluded ({reason})")
                    log_file.write(f"{record.id}\tExcluded\t{reason}\t{hit_def} ~{int(identity_pct)}%\n")
            else:
                reason = "No alignment"
                print(f"  ❌ {reason}")
                log_file.write(f"{record.id}\tExcluded\t{reason}\tN/A\n")

            processed_ids.add(record.id)
            time.sleep(3)

    # Write Summary CSV
    with open(summary_csv_path, "w", newline='', encoding='utf-8') as csvfile:
        fieldnames = ["Species", "Count"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for species, count in sorted(species_counts.items(), key=lambda x: x[1], reverse=True):
            writer.writerow({"Species": species, "Count": count})

    print(f"\n📊 Summary CSV saved: {summary_csv_path}")
    print(f"✅ Process for barcode{barcode} completed.\n")

    return True

def main():
    start_barcode_input = input("Enter the barcode number (e.g., 05 for barcode05): ").strip().zfill(2)
    try:
        start_barcode_num = int(start_barcode_input)
    except ValueError:
        sys.exit("❌ Invalid barcode number entered.")

    barcode_folders = get_barcode_folders(base_dir)
    if not barcode_folders:
        sys.exit(f"❌ No barcode folders found in {base_dir}")

    # Filter folders to those with barcode >= start_barcode_num
    folders_to_process = [(num, path) for num, path in barcode_folders if num >= start_barcode_num]

    if not folders_to_process:
        sys.exit(f"❌ No barcode folders found starting from barcode{start_barcode_input}")

    print(f"\nFound barcode folders to process: {[f'barcode{str(n).zfill(2)}' for n, _ in folders_to_process]}")

    for num, folder_path in folders_to_process:
        run_blast_for_barcode(num, folder_path)

    print("🚀 All done!")

if __name__ == "__main__":
    main()
