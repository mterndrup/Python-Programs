#!/usr/bin/env python3
"""
Recursively extract .gz, .tar.gz/.tgz and .zip files
inside the Round2 directory, leaving the results
in the folder where each archive is located.

Author: you 😊
"""

from pathlib import Path
import gzip
import shutil
import tarfile
import zipfile

# -------- EDIT THIS LINE IF YOUR PATH CHANGES ----------
BASE_DIR = Path(r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\Round2")
# -------------------------------------------------------

def extract_gzip(gz_path: Path):
    """Decompress a single-file .gz archive (e.g. *.fastq.gz)."""
    out_path = gz_path.with_suffix('')  # remove the .gz
    if out_path.exists():
        print(f"  • Skipping (already exists): {out_path.name}")
        return
    print(f"  • Extracting {gz_path.name} → {out_path.name}")
    with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def extract_tar_gz(tar_path: Path):
    """Extract .tar.gz or .tgz archives."""
    print(f"  • Extracting TAR-GZ: {tar_path.name}")
    with tarfile.open(tar_path, 'r:gz') as tar:
        tar.extractall(tar_path.parent)

def extract_zip(zip_path: Path):
    """Extract .zip archives."""
    print(f"  • Extracting ZIP: {zip_path.name}")
    with zipfile.ZipFile(zip_path, 'r') as z:
        z.extractall(zip_path.parent)

def is_tar_gz(p: Path) -> bool:
    """True if file looks like *.tar.gz or *.tgz."""
    return (p.suffix.lower() == ".tgz") or (p.suffixes[-2:] == ['.tar', '.gz'])

def main():
    if not BASE_DIR.exists():
        print(f"Base directory not found: {BASE_DIR}")
        return

    print(f"Scanning {BASE_DIR}…")

    for path in BASE_DIR.rglob("*"):
        if path.is_file():
            lower = path.name.lower()
            try:
                if lower.endswith(".zip"):
                    extract_zip(path)
                elif is_tar_gz(path):
                    extract_tar_gz(path)
                elif lower.endswith(".gz"):
                    extract_gzip(path)
            except Exception as exc:
                print(f"  ! Error processing {path}: {exc}")

    print("Done!")

if __name__ == "__main__":
    main()
