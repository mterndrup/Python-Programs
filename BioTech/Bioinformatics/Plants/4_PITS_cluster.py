import subprocess
import os

# === Set vsearch path ===
vsearch_path = (
    r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
    r"\Bioinformatics\Tools\VSEARCH\vsearch-2.30.0-win-x86_64\bin\vsearch.exe"
)

def process_barcode(base_path, barcode_num, region):
    """Process a single barcode for the given region (ITS2 or LSU)."""
    # Construct input paths
    input_fasta = os.path.join(
        base_path,
        f"barcode{barcode_num}",
        "ITSx-All",
        f"barcode{barcode_num}_combined_ITSx.{region}.fasta"
    )

    if not os.path.isfile(input_fasta):
        print(f"Error: FASTA file not found:\n{input_fasta}")
        return

    # Setup output paths
    base_dir = os.path.dirname(input_fasta)
    base_name = os.path.splitext(os.path.basename(input_fasta))[0]

    derep_fasta = os.path.join(base_dir, f"{base_name}_derep.fasta")
    centroids_fasta = os.path.join(base_dir, f"{base_name}_centroids.fasta")
    uc_file = os.path.join(base_dir, f"{base_name}_clusters.uc")
    dominant_fasta = os.path.join(base_dir, f"{base_name}_dominant.fasta")

    # Step 1: Dereplicate
    print(f"\nProcessing barcode {barcode_num} ({region}): Step 1 - Dereplicating sequences...")
    subprocess.run([
        vsearch_path,
        "--derep_fulllength", input_fasta,
        "--output", derep_fasta,
        "--sizeout"
    ], check=True)

    # Step 2: Cluster at 99%
    print(f"Processing barcode {barcode_num} ({region}): Step 2 - Clustering sequences at 99% identity...")
    subprocess.run([
        vsearch_path,
        "--cluster_fast", derep_fasta,
        "--id", "0.99",
        "--strand", "both",
        "--centroids", centroids_fasta,
        "--uc", uc_file,
        "--sizeout"
    ], check=True)

    # Step 3: Extract top (dominant) sequence
    print(f"Processing barcode {barcode_num} ({region}): Step 3 - Extracting the most abundant sequence...")
    top_seq = ""
    top_seq_id = ""
    with open(centroids_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                top_seq_id = line.strip()
                continue
            top_seq = line.strip()
            break

    with open(dominant_fasta, "w") as out:
        out.write(f"{top_seq_id}\n{top_seq}\n")

    print(f"Dominant sequence written to {dominant_fasta}")


def main():
    base_path = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus2"

    # Ask whether to process one barcode or all
    choice = input("Process one barcode or all? Enter 'one' or 'all': ").strip().lower()
    if choice not in ["one", "all"]:
        print("Invalid choice. Please enter 'one' or 'all'.")
        return

    # Ask for ITS2 or LSU
    region = input("Which region to process? Enter 'ITS2' or 'LSU': ").strip().upper()
    if region not in ["ITS2", "LSU"]:
        print("Invalid choice. Must be 'ITS2' or 'LSU'.")
        return

    if choice == "one":
        barcode_num = input("Enter barcode number (e.g., 73): ").strip()
        if not barcode_num.isdigit():
            print("Invalid barcode number. Must be digits only.")
            return
        barcode_num = barcode_num.zfill(2)
        process_barcode(base_path, barcode_num, region)

    elif choice == "all":
        all_barcodes = [d for d in os.listdir(base_path) if d.startswith("barcode")]
        if not all_barcodes:
            print("No barcode directories found.")
            return
        print(f"Found {len(all_barcodes)} barcodes. Processing all ({region})...")
        for barcode_dir in sorted(all_barcodes):
            barcode_num = barcode_dir.replace("barcode", "")
            process_barcode(base_path, barcode_num, region)


if __name__ == "__main__":
    main()
