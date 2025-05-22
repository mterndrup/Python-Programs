import subprocess
import os

def run_mafft():
    # Hardcoded input/output paths (adjusted for WSL)
    input_fasta = ("/mnt/c/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/"
                   "BioTech/Bioinformatics/Chaparral/wildfirePlants-DNA-nanopore-sequence/"
                   "fastq_pass/barcode47/barcode47_PITS_combine_cluster.fasta")
    output_fasta = ("/mnt/c/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/"
                    "Bioinformatics/Chaparral/wildfirePlants-DNA-nanopore-sequence/"
                    "fastq_pass/barcode47/barcode47_PITS_combine_cluster_aligned.fasta")
    threads = 4

    if not os.path.isfile(input_fasta):
        print(f"Error: Input file '{input_fasta}' not found.")
        return

    cmd = [
        "mafft",
        "--thread", str(threads),
        "--auto",
        input_fasta
    ]

    print("Running MAFFT command:")
    print(" ".join(cmd))

    try:
        with open(output_fasta, "w") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            print(f"Alignment completed successfully! Output saved to {output_fasta}")
        else:
            print("MAFFT error:")
            print(result.stderr)
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    run_mafft()