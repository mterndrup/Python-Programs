import subprocess

def run_dada2(base_path, r_script_path, rscript_exe):
    cmd = [
        rscript_exe,
        r_script_path,
        base_path
    ]

    print(f"Running command:\n{' '.join(cmd)}")

    try:
        # Capture stdout and stderr for debugging
        result = subprocess.run(cmd, check=True, text=True, capture_output=False)
        print("=== STDOUT ===")
        print(result.stdout)
        print("=== STDERR ===")
        print(result.stderr)
        print("DADA2 script completed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error: DADA2 script failed.")
        print("Return code:", e.returncode)
        print("=== STDOUT ===")
        print(e.stdout)
        print("=== STDERR ===")
        print(e.stderr)

if __name__ == "__main__":
    base_path = r"/BioTech/Bioinformatics/Chaparral/Round2"
    r_script_path = r"/BioTech/Bioinformatics/Tools/DADA2/run_dada2.R"
    rscript_exe = r"C:\Program Files\R\R-4.5.0\bin\x64\Rscript.exe"

    run_dada2(base_path, r_script_path, rscript_exe)
