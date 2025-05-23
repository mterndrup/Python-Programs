from collections import Counter
import re


def summarize_fasta_headers_ignore_percentage():
    fasta_file = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                  r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
                  r"\fastq_pass\barcode52\barcode52_FITS_fungi_only.fasta")
    output_file = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                   r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
                   r"\fastq_pass\barcode52\summary_counts.txt")

    counts = Counter()

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                header_clean = re.sub(r'\s*~.*', '', header)
                counts[header_clean] += 1

    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)

    with open(output_file, 'w') as out:
        for header, count in sorted_counts:
            line = f"{header}: {count}"
            print(line)  # print to terminal
            out.write(line + "\n")  # write to file


if __name__ == "__main__":
    summarize_fasta_headers_ignore_percentage()
