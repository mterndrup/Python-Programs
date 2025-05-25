import subprocess
import os

def add_barcode_prefix_to_fasta(fasta_file, barcode_prefix):
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    with open(fasta_file, 'w') as f:
        for line in lines:
            if line.startswith('>'):
                line = f">{barcode_prefix}_{line[1:]}"
            f.write(line)

def add_barcode_prefix_to_uc(uc_file, barcode_prefix):
    with open(uc_file, 'r') as f:
        lines = f.readlines()

    with open(uc_file, 'w') as f:
        for line in lines:
            if line.startswith('#') or line.strip() == '':
                # keep comments or empty lines as is
                f.write(line)
                continue
            parts = line.strip().split('\t')
            # Add prefix to the ID in 9th column if it exists
            if len(parts) > 8:
                parts[8] = f"{barcode_prefix}_{parts[8]}"
            f.write('\t'.join(parts) + '\n')

def cluster_its_reads(input_fasta, identity=0.99, threads=4):
    base = os.path.splitext(input_fasta)[0]
    centroids_out = f"{base}_cluster.fasta"
    uc_out = f"{base}_cluster.uc"

    cmd = [
        "vsearch",
        "--cluster_fast", input_fasta,
        "--id", str(identity),
        "--centroids", centroids_out,
        "--uc", uc_out,
        "--threads", str(threads)
    ]

    print("Running command:")
    print(" ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print("Clustering complete!")
        print(f"Centroids saved to: {centroids_out}")
        print(f"Cluster map saved to: {uc_out}")

        barcode = None
        for part in input_fasta.split(os.sep):
            if part.startswith("barcode"):
                barcode = part
                break
        if not barcode:
            filename = os.path.basename(input_fasta)
            if "barcode" in filename:
                start = filename.index("barcode")
                barcode = filename[start:start+8]
            else:
                barcode = "unknownbarcode"

        print(f"Adding barcode prefix '{barcode}' to cluster fasta headers...")
        add_barcode_prefix_to_fasta(centroids_out, barcode)
        print("Header prefixing complete.")

        print(f"Adding barcode prefix '{barcode}' to .uc file sequence IDs...")
        add_barcode_prefix_to_uc(uc_out, barcode)
        print(".uc file prefixing complete.")

    else:
        print("Error during clustering:")
        print(result.stderr)


if __name__ == "__main__":
    input_file = ("/mnt/c/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/"
                  "BioTech/Bioinformatics/Chaparral/wildfirePlants-DNA-nanopore-"
                  "sequence/fastq_pass/barcode48/barcode48_PITS_combine.fasta")

    if not os.path.isfile(input_file):
        print(f"Error: File '{input_file}' not found.")
        exit(1)

    cluster_its_reads(input_file)