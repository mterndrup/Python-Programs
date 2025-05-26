import os
import subprocess
from Bio import SeqIO

# === User parameters ===
main_folder = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
               r"\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence"
               r"\fastq_pass")

vsearch_path = (r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech"
                r"\Bioinformatics\Chaparral\vsearch-2.30.0-win-x86_64"
                r"\vsearch-2.30.0-win-x86_64\bin\vsearch.exe")


# === Helper function to sort FASTA by sequence length descending ===
def sort_fasta_by_length(input_fasta, sorted_fasta):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    records.sort(key=lambda r: len(r.seq), reverse=True)  # longest first
    SeqIO.write(records, sorted_fasta, "fasta")
    print(f"🔀 Sorted input FASTA by sequence length and wrote to: {sorted_fasta}")


# === Function to process one barcode folder ===
def process_barcode_folder(barcode_folder):
    print(f"\n➡️ Processing folder: {barcode_folder}")

    # Find the FASTA file in this folder (assuming one .fasta file)
    fasta_files = [f for f in os.listdir(barcode_folder) if f.endswith(".fasta")]
    if not fasta_files:
        print("⚠️ No FASTA file found, skipping.")
        return

    input_fasta = os.path.join(barcode_folder, fasta_files[0])

    # Create VSEARCH output folder inside barcode_folder
    vsearch_folder = os.path.join(barcode_folder, "VSEARCH")
    os.makedirs(vsearch_folder, exist_ok=True)

    # Create a sorted fasta file path inside VSEARCH folder
    sorted_fasta = os.path.join(vsearch_folder, "sorted_" + fasta_files[0])

    # Sort the input fasta by read length descending
    sort_fasta_by_length(input_fasta, sorted_fasta)

    prefix_name = os.path.splitext(fasta_files[0])[0]
    cluster_output_prefix = os.path.join(vsearch_folder, prefix_name + "_clusters")

    # Step 1: Run vsearch clustering on the sorted fasta file
    print("🔧 Running vsearch clustering...")

    cluster_cmd = [
        vsearch_path,
        "--cluster_fast", sorted_fasta,
        "--id", "0.97",  # keep 97% identity for species-level resolution
        "--strand", "both",  # consider both strands
        "--sizeout",  # add cluster sizes to fasta headers
        "--centroids", cluster_output_prefix + "_centroids.fasta",
        "--uc", cluster_output_prefix + ".uc",
        "--threads", "4"
    ]

    try:
        subprocess.run(cluster_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ vsearch clustering failed: {e}")
        return

    print("✅ Clustering finished.")

    # Step 2: Parse .uc file to build cluster map
    uc_file = cluster_output_prefix + ".uc"
    cluster_map = {}

    print(f"📖 Parsing cluster file: {uc_file}")

    with open(uc_file, "r") as uc:
        for line in uc:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            record_type = parts[0]
            query_id = parts[8]

            if record_type == "S":
                cluster_map[query_id] = [query_id]
            elif record_type == "H":
                centroid_id = parts[9]
                if centroid_id not in cluster_map:
                    cluster_map[centroid_id] = [centroid_id]
                cluster_map[centroid_id].append(query_id)

    print(f"🔎 Parsed {len(cluster_map)} clusters.")

    # Step 3: Load sequences from sorted FASTA
    print("📥 Loading input FASTA sequences...")
    seq_dict = SeqIO.to_dict(SeqIO.parse(sorted_fasta, "fasta"))
    print(f"Loaded {len(seq_dict)} sequences.")

    # Step 4: Write clustered FASTA with combined headers inside VSEARCH folder
    output_clustered_fasta = cluster_output_prefix + "_clustered.fasta"
    print(f"📁 Writing clustered FASTA to: {output_clustered_fasta}")

    with open(output_clustered_fasta, "w") as out_fa:
        for centroid_id, members in cluster_map.items():
            if len(members) == 1:
                combined_header = centroid_id
            else:
                combined_header = f"{centroid_id}|{'/'.join([m for m in members if m != centroid_id])}"

            if centroid_id in seq_dict:
                seq_record = seq_dict[centroid_id]
                seq_record.id = combined_header
                seq_record.description = ""
                SeqIO.write(seq_record, out_fa, "fasta")
            else:
                print(f"⚠️ Warning: centroid sequence {centroid_id} not found in FASTA!")

    print("🎉 Done with this folder.")


# === Main loop through barcode folders ===
for folder_name in os.listdir(main_folder):
    folder_path = os.path.join(main_folder, folder_name)
    if os.path.isdir(folder_path) and folder_name.startswith("barcode"):
        process_barcode_folder(folder_path)

print("\n🏁 All barcode folders processed.")
