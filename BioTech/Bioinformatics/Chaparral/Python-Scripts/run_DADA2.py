import subprocess
import os
from datetime import datetime

def run_dada2(base_path, r_script_path, rscript_exe, multithread=True):
    """
    Run DADA2 R script for _cutadapt_combined.fasta files with robust logging and recursive detection.
    """
    # Check paths
    if not os.path.exists(base_path):
        raise FileNotFoundError(f"Base path does not exist: {base_path}")
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"R script not found: {r_script_path}")
    if not os.path.exists(rscript_exe):
        raise FileNotFoundError(f"Rscript executable not found: {rscript_exe}")

    # Detect barcode folders
    barcode_dirs = [
        os.path.join(base_path, d)
        for d in os.listdir(base_path)
        if os.path.isdir(os.path.join(base_path, d)) and "barcode" in d.lower()
    ]

    # Filter barcode directories with numeric > 19
    filtered_barcodes = []
    for d in barcode_dirs:
        name = os.path.basename(d)
        num_part = ''.join(filter(str.isdigit, name))
        if num_part and int(num_part) > 19:
            filtered_barcodes.append(d)

    if not filtered_barcodes:
        print("No barcode directories found with numeric > 19.")
        return

    print("Filtered barcode directories:", [os.path.basename(d) for d in filtered_barcodes])

    # Check for _cutadapt_combined.fasta files recursively
    barcodes_with_fasta = []
    for d in filtered_barcodes:
        fasta_files = []
        for root, _, files in os.walk(d):
            for f in files:
                if f.endswith("_cutadapt_combined.fasta"):
                    fasta_files.append(os.path.join(root, f))
        if fasta_files:
            barcodes_with_fasta.append(d)
            print(f"{os.path.basename(d)}: Found FASTA file {os.path.basename(fasta_files[0])} in {os.path.relpath(os.path.dirname(fasta_files[0]), d)}")
        else:
            print(f"{os.path.basename(d)}: No _cutadapt_combined.fasta found, skipping.")

    if not barcodes_with_fasta:
        print("No FASTA files found in any barcode directories.")
        return

    # Prepare Rscript command
    cmd = [rscript_exe, r_script_path, base_path]
    if not multithread:
        cmd.append("singlethread")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(base_path, f"dada2_wrapper_log_{timestamp}.txt")

    print("\nRunning R script...")
    print("Command:", ' '.join(cmd))
    print(f"Logging output to: {log_file}\n")

    with open(log_file, "w", encoding="utf-8") as log:
        try:
            result = subprocess.run(
                cmd,
                text=True,
                capture_output=True
            )

            # Write stdout/stderr to log
            log.write("=== STDOUT ===\n")
            log.write(result.stdout or "None\n")
            log.write("\n=== STDERR ===\n")
            log.write(result.stderr or "None\n")
            log.write(f"\nReturn code: {result.returncode}\n")

            # Print stdout/stderr to console
            print("=== STDOUT ===")
            print(result.stdout or "None")
            print("=== STDERR ===")
            print(result.stderr or "None")
            print(f"Return code: {result.returncode}")

            if result.returncode == 0:
                print("\nDADA2 script completed successfully.")
            else:
                print("\nDADA2 script failed. Check log file for details.")

        except Exception as e:
            log.write(f"Exception occurred: {str(e)}\n")
            print(f"Exception occurred: {str(e)}")

if __name__ == "__main__":
    # Update these paths to your environment
    base_path = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus2"  # parent folder of barcode dirs
    r_script_path = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Tools\DADA2\run_dada2.R"
    rscript_exe = r"C:\Program Files\R\R-4.5.0\bin\x64\Rscript.exe"

    # Run DADA2
    run_dada2(base_path, r_script_path, rscript_exe, multithread=False)
