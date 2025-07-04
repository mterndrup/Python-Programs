from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import time
import csv
import os
import re

Entrez.email = "your_email@example.com"  # Replace with your email
search_term = "Banisteriopsis caapi[Organism]"
max_results = 1000

output_folder = "References"
os.makedirs(output_folder, exist_ok=True)


def extract_region_status(description, region_name):
    # Check standard partial/complete mention
    match = re.search(rf"{region_name}.*?(partial|complete)", description, re.IGNORECASE)
    if match:
        return match.group(1).capitalize()

    # NEW: If it's in a "genes for ..." block, mark as "Present"
    genes_for_block = re.search(r"genes for (.+?)\s*(?:partial|complete|sequence|$)", description, re.IGNORECASE)
    if genes_for_block:
        gene_list = genes_for_block.group(1).lower()
        if region_name.lower() in gene_list:
            return "Present"

    # Standard fallback
    elif region_name.lower() in description.lower():
        return "Present"

    return "-"


def extract_function(record, its1, r5_8s, its2, lsu, r18s, r26s):
    # Check if it's a complete genome
    if "complete genome" in record.description.lower():
        return "Complete Genome"

    # Check for known rRNA regions directly in the description text
    description_text = record.description.lower()
    rRNA_keywords = ["18s", "its1", "5.8s", "its2", "26s"]

    if any(keyword in description_text for keyword in rRNA_keywords):
        return "Identifier"

    # Fallback: look for annotations
    for feature in record.features:
        for key in ['product', 'note', 'gene']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]

    return record.description


def has_protein(record):
    # Check if there is a protein-coding sequence (CDS) with translation
    for feature in record.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            return True
    return False


def save_fasta(record, filename):
    # Only save if the file doesn't already exist
    if not os.path.exists(filename):
        with open(filename, "w") as f:
            SeqIO.write(record, f, "fasta")


def fetch_accession_ids(term, max_results):
    handle = Entrez.esearch(db="nucleotide", term=term, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]


def fetch_metadata(ncbi_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        accession = record.id
        description = record.description
        strain = "-"
        for feature in record.features:
            if feature.type == "source" and "strain" in feature.qualifiers:
                strain = feature.qualifiers["strain"][0]
                break

        # Extract the ribosomal RNA region status
        its1 = extract_region_status(description, "internal transcribed spacer 1")
        r5_8s = extract_region_status(description, "5.8S ribosomal RNA")
        its2 = extract_region_status(description, "internal transcribed spacer 2")
        lsu = extract_region_status(description, "large subunit ribosomal RNA")
        r18s = extract_region_status(description, "18S ribosomal RNA")  # Check for 18S
        r26s = extract_region_status(description, "26S ribosomal RNA")  # Check for 26S

        # Determine the function based on the regions and other features
        function = extract_function(record, its1, r5_8s, its2, lsu, r18s, r26s)
        protein_flag = "Yes" if has_protein(record) else "-"

        # Save FASTA only if it doesn't already exist
        fasta_filename = os.path.join(output_folder, f"{accession}.fasta")
        save_fasta(record, fasta_filename)

        # Return all the information, including the region statuses
        return accession, strain, its1, r5_8s, its2, lsu, r18s, r26s, function, protein_flag, description

    except Exception as e:
        return "Error", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", str(e)


# Run the script to fetch metadata and save results
ids = fetch_accession_ids(search_term, max_results)
results = []

for idx, ncbi_id in enumerate(ids):
    row = fetch_metadata(ncbi_id)
    print(f"{idx + 1}. {row[0]} | Protein: {row[9]} | Function: {row[8]}")
    results.append(row)
    time.sleep(0.4)

# Save to CSV with columns for the ribosomal RNA regions and function
with open("Aya/References/b.caapi-genetics.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Accession", "Strain",
        "ITS1", "5.8S", "ITS2", "LSU", "18S", "26S",  # Only include 26S here
        "Function", "Protein", "Full Description"
    ])
    writer.writerows(results)