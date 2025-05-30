from collections import Counter
import re
import csv

def summarize_fasta_headers_ignore_percentage():
    fasta_file = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                  r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
                  r"\fastq_pass\barcode08\barcode08_FITS_fungi_only.fasta")
    output_file = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                   r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
                   r"\fastq_pass\barcode08\barcode08_sum.csv")  # keep .csv extension

    counts = Counter()

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                header_clean = re.sub(r'\s*~.*', '', header)
                counts[header_clean] += 1

    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Species', 'Count'])  # header row with your desired names
        for header, count in sorted_counts:
            writer.writerow([header, count])
            print(f"{header}: {count}")  # print to terminal

if __name__ == "__main__":
    summarize_fasta_headers_ignore_percentage()
