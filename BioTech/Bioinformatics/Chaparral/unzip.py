from pathlib import Path
import gzip
import shutil
import tarfile
import zipfile
from datetime import datetime

BASE_DIR = Path(r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Wild")
LOG_DIR = Path(r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Logs")

def log_message(log_file, message):
    print(message)
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(message + "\n")

def extract_gzip(gz_path: Path, raw_dir: Path, log_file: Path):
    out_path = raw_dir / gz_path.with_suffix('').name  # Remove .gz suffix
    if out_path.exists():
        log_message(log_file, f"  • Skipping (exists): {out_path.name}")
        return
    log_message(log_file, f"  • Extracting {gz_path.name} → {out_path.name}")
    with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    shutil.move(str(gz_path), str(raw_dir / gz_path.name))

def extract_tar_gz(tar_path: Path, raw_dir: Path, log_file: Path):
    log_message(log_file, f"  • Extracting TAR-GZ: {tar_path.name}")
    with tarfile.open(tar_path, 'r:gz') as tar:
        tar.extractall(raw_dir)
    shutil.move(str(tar_path), str(raw_dir / tar_path.name))

def extract_zip(zip_path: Path, raw_dir: Path, log_file: Path):
    log_message(log_file, f"  • Extracting ZIP: {zip_path.name}")
    with zipfile.ZipFile(zip_path, 'r') as z:
        z.extractall(raw_dir)
    shutil.move(str(zip_path), str(raw_dir / zip_path.name))

def is_tar_gz(p: Path) -> bool:
    return (p.suffix.lower() == ".tgz") or (p.suffixes[-2:] == ['.tar', '.gz'])

def main():
    if not BASE_DIR.exists():
        print(f"Base directory not found: {BASE_DIR}")
        return

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = LOG_DIR / f"extraction_log_{timestamp}.txt"

    start_time = datetime.now()
    log_message(log_file, f"Script started at {start_time}")
    log_message(log_file, f"Scanning base directory: {BASE_DIR}")

    total_processed = 0
    processed_gz_files = set()  # Track extracted gz paths

    # Iterate ONLY immediate subfolders (barcode folders)
    for barcode_dir in BASE_DIR.iterdir():
        if not barcode_dir.is_dir():
            continue

        # Skip any folder named unclassified (case insensitive)
        if barcode_dir.name.lower() == "unclassified":
            log_message(log_file, f"Skipping 'unclassified' folder: {barcode_dir}")
            continue

        raw_dir = barcode_dir / "Raw"
        raw_dir.mkdir(exist_ok=True)

        # Process files ONLY in barcode_dir (exclude Raw folder)
        for file_path in barcode_dir.iterdir():
            if file_path.is_file():
                lower_name = file_path.name.lower()

                # Skip files inside Raw or any other folder (should be none here)
                # We only scan files directly inside barcode_dir

                try:
                    if lower_name.endswith(".zip"):
                        extract_zip(file_path, raw_dir, log_file)
                        total_processed += 1
                    elif is_tar_gz(file_path):
                        extract_tar_gz(file_path, raw_dir, log_file)
                        total_processed += 1
                    elif lower_name.endswith(".gz"):
                        # Don't extract again if already done
                        if file_path in processed_gz_files:
                            log_message(log_file, f"  • Already processed {file_path.name}, skipping.")
                            continue
                        if not is_tar_gz(file_path):
                            extract_gzip(file_path, raw_dir, log_file)
                            total_processed += 1
                            processed_gz_files.add(file_path)
                except Exception as e:
                    log_message(log_file, f"  ! Error processing {file_path}: {e}")

    log_message(log_file, f"Done! Total files processed: {total_processed}")

    end_time = datetime.now()
    duration = end_time - start_time
    log_message(log_file, f"Script finished at {end_time}")
    log_message(log_file, f"Duration: {duration}")

if __name__ == "__main__":
    main()
