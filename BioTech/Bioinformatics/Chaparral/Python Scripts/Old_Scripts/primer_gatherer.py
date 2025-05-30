from Bio import SeqIO
from Bio.Seq import Seq
import os
import re

# Define ITS primers with names and orientations
ITS_PRIMERS = {
    "ITS-S2F (forward)": "ATGCGATACTTGGTGTGAAT",
    "ITS4 (forward)": "TCCTCCGCTTATTGATATGC",
    "ITS-S2F (reverse complement)": "ATTCACACCAAGTATCGCAT",
    "ITS4 (reverse complement)": "GCATATCAATAAGCGGAGGA"
}

def find_primer(sequence, primers_dict):
    for name, primer_seq in primers_dict.items():
        if primer_seq in sequence:
            return name, primer_seq
    return None, None

# Set working directory
working_dir = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs"
               r"\BioTech\Bioinformatics\Tribal-Medicine\Cacti\DNA_Sequences\barcode48")

# Get folder name from working_dir path
folder_name = os.path.basename(working_dir.rstrip(os.sep))

output_fasta = os.path.join(working_dir, f"{folder_name}_PITS_combine.fasta")

total_count = 0
with open(output_fasta, "w") as out_handle:
    for filename in os.listdir(working_dir):
        if filename.endswith("_trimmed.fastq"):
            # Extract the numeric trimmed label (e.g. 0_trimmed) from filename
            match = re.search(r"(\d+_trimmed)\.fastq$", filename)
            if match:
                trimmed_label = match.group(1)
            else:
                trimmed_label = os.path.splitext(filename)[0]

            input_fastq = os.path.join(working_dir, filename)
            count = 0
            print(f"Processing {filename}...")
            for record in SeqIO.parse(input_fastq, "fastq"):
                primer_name, primer_seq = find_primer(str(record.seq), ITS_PRIMERS)
                if primer_name:
                    print(f"🔍 {trimmed_label}")
                    print(f"{primer_name}: {primer_seq}")
                    SeqIO.write(record, out_handle, "fasta")
                    count += 1
            print(f"✅ {trimmed_label}: {count} reads added.")
            total_count += count

print(f"\n🎉 All done! Total {total_count} reads written to:\n{output_fasta}")