import subprocess
import os

def windows_to_wsl_path(win_path):
    drive, path_rest = os.path.splitdrive(win_path)
    drive_letter = drive.rstrip(':').lower()
    # Remove colon from drive and ensure forward slashes
    path_rest = path_rest.replace('\\', '/')
    # Remove leading slash if present in path_rest to avoid double slashes
    if path_rest.startswith('/'):
        path_rest = path_rest[1:]
    return f"/mnt/{drive_letter}/{path_rest}"

def run_nanoclust_wsl(fasta_input_win, output_dir_win, conda_env="nanoclust", cpu=4):
    fasta_input_wsl = windows_to_wsl_path(fasta_input_win)
    output_dir_wsl = windows_to_wsl_path(output_dir_win)

    cmd = (
        f'eval "$(conda shell.bash hook)" && '
        f'conda activate {conda_env} && '
        f'NanoCLUST -i "{fasta_input_wsl}" -o "{output_dir_wsl}" --cpu {cpu}'
    )

    print("🔹 Running NanoCLUST inside WSL environment with command:")
    print(cmd)

    subprocess.run(["bash", "-c", cmd], check=True)

def main():
    fasta_input = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence\fastq_pass\barcode05\barcode05_FITS_full.fasta"
    output_dir = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence\fastq_pass\barcode05\NanoCLUST_results"

    run_nanoclust_wsl(fasta_input, output_dir, conda_env="nanoclust", cpu=4)

if __name__ == "__main__":
    main()
