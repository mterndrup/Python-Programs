import os
import subprocess

# Change this to your parent directory containing barcode folders
parent_dir = "/mnt/c/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/Bioinformatics/Chaparral/wildfirePlants-DNA-nanopore-sequence/fastq_pass"

# Path to ITSx executable script
itsx_path = "/mnt/c/Users/ketgl/Downloads/ITSx_1.1.3/ITSx_1.1.3/ITSx"

# Number of CPUs to use
cpus = "4"

for foldername in os.listdir(parent_dir):
    folderpath = os.path.join(parent_dir, foldername)
    if os.path.isdir(folderpath) and foldername.startswith("barcode"):
        print(f"Processing folder: {foldername}")

        itsx_output_dir = os.path.join(folderpath, "ITSx")
        os.makedirs(itsx_output_dir, exist_ok=True)

        # Find fasta files matching pattern *_FITS_final_trimmed.fasta
        for filename in os.listdir(folderpath):
            if filename.endswith("_combined.fasta"):
                input_fasta = os.path.join(folderpath, filename)
                output_prefix = os.path.join(itsx_output_dir, filename.replace(".fasta", "_ITSx"))

                cmd = [
                    "perl",
                    itsx_path,
                    "-i", input_fasta,
                    "-o", output_prefix,
                    "--cpu", cpus,
                    "--preserve", "T",
                    "--allow_partial",
                    "--save_regions", "all"
                ]

                print("Running ITSx for:", input_fasta)
                try:
                    subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running ITSx on {input_fasta}: {e}")
