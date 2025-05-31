import os
import re

cds_file = 'OR574827.1-all-orfs-DNA.cds'
output_folder = 'fastas'

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def write_fasta(record_id, sequence, folder):
    fasta_path = os.path.join(folder, f"{record_id}.fasta")
    with open(fasta_path, 'w') as f:
        f.write(f">{record_id}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i + 60] + '\n')


def write_gff3(record_id, accession, start, end, strand, orf_id, folder):
    gff_path = os.path.join(folder, f"{record_id}.gff3")
    with open(gff_path, 'w') as gff:
        gff.write("##gff-version 3\n")
        gff.write(f"{accession}\tcds_parser\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={orf_id}\n")
        gff.write(f"{accession}\tcds_parser\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={orf_id}_cds;Parent={orf_id}\n")


with open(cds_file, 'r') as f:
    content = f.read()

records = content.strip().split('>')[1:]  # Skip the initial empty split
for record in records:
    record = record.strip()
    match = re.match(r"gi\|\d+\|\w+\|([\w.]+)\|\:([c]?\d+)-(\d+)\s+(ORF\d+)\s+CDS\n(.+)", record, re.DOTALL)

    if match:
        accession = match.group(1)
        start = match.group(2)
        end = match.group(3)
        orf_id = match.group(4)
        sequence = match.group(5).replace('\n', '').replace(' ', '')

        # Determine strand
        strand = '-' if start.startswith('c') else '+'
        start = start.lstrip('c')

        record_id = f"{accession}_{orf_id}"
        accession_folder = os.path.join(output_folder, accession)
        if not os.path.exists(accession_folder):
            os.makedirs(accession_folder)

        # Write the outputs
        write_fasta(record_id, sequence, accession_folder)
        write_gff3(record_id, accession, start, end, strand, orf_id, accession_folder)

    else:
        print(f"⚠️ Skipping invalid record (no match): {record[:60]}...")