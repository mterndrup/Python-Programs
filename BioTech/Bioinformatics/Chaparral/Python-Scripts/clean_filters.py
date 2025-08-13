import os

def fastq_to_fasta_keep_headers(fastq_path, fasta_path):
    with open(fastq_path, 'r') as fq, open(fasta_path, 'w') as fa:
        while True:
            header = fq.readline().rstrip()
            if not header:
                break
            seq = fq.readline().rstrip()
            plus = fq.readline().rstrip()
            qual = fq.readline().rstrip()

            if not (header.startswith('@') and plus.startswith('+')):
                raise ValueError("File not in proper FASTQ format")

            fa.write('>' + header[1:] + '\n')
            fa.write(seq + '\n')

base_folder = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\Round2"

# Loop through each item in the base folder
for entry in os.listdir(base_folder):
    subfolder_path = os.path.join(base_folder, entry)
    # Check if it is a directory (barcode folder)
    if os.path.isdir(subfolder_path):
        # List all files in the barcode folder
        for filename in os.listdir(subfolder_path):
            if filename.endswith(".fastq") or filename.endswith(".fq"):
                fastq_file = os.path.join(subfolder_path, filename)
                fasta_file = os.path.join(subfolder_path, filename.rsplit('.', 1)[0] + ".fasta")
                fastq_to_fasta_keep_headers(fastq_file, fasta_file)
                print(f"Converted {fastq_file} to {fasta_file}")
