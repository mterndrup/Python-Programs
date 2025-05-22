import os
import subprocess
import logging
from datetime import datetime

# Setup logging
log_file = os.path.expanduser("~/scripts/porechop_run.log")
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

def main():
    base_dir = "/mnt/c/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/Bioinformatics/Chaparral/wildfirePlants-DNA-nanopore-sequence/fastq_pass"

    print("🔪 Starting Porechop batch processing")
    logging.info("Starting Porechop batch processing")

    for barcode_folder in sorted(os.listdir(base_dir)):
        barcode_path = os.path.join(base_dir, barcode_folder)
        if not os.path.isdir(barcode_path) or not barcode_folder.startswith("barcode"):
            continue

        print(f"📂 Processing folder: {barcode_folder}")
        logging.info(f"Processing folder: {barcode_folder}")

        for file in os.listdir(barcode_path):
            if file.endswith(".fastq"):
                input_file = os.path.join(barcode_path, file)
                filename_no_ext, ext = os.path.splitext(file)
                output_file = os.path.join(barcode_path, f"{filename_no_ext}_trimmed{ext}")

                print(f"  🔪 Running Porechop on {file}...")
                logging.info(f"Running Porechop on {file}")

                cmd = [
                    "porechop",
                    "-i", input_file,
                    "-o", output_file,
                    "--threads", "4"
                ]

                try:
                    subprocess.run(cmd, check=True)
                    print(f"  ✅ Successfully processed {file}")
                    logging.info(f"Successfully processed {file}")
                except subprocess.CalledProcessError as e:
                    print(f"  ❌ Error running Porechop on {file}: {e}")
                    logging.error(f"Error running Porechop on {file}: {e}")

    print("🔪 Completed Porechop batch processing")
    logging.info("Completed Porechop batch processing")

if __name__ == "__main__":
    main()