from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import time
import csv
import os
import re

Entrez.email = "your_email@example.com"  # Replace with your email
search_term = "Pestalotiopsis microspora[Organism]"
max_results = 1000

output_folder = "DNA-Sequences"
os.makedirs(output_folder, exist_ok=True)

def extract_region_status(description, region_name):
    match = re.search(rf"{region_name}.*?(partial|complete)", description, re.IGNORECASE)
    if match:
        return match.group(1).capitalize()
    elif region_name.lower() in description.lower():
        return "Present"
    else:
        return "-"

def extract_function(record, its1, r5_8s, its2, lsu):
    if any(region != "-" for region in [its1, r5_8s, its2, lsu]):
        return "Identifier"
    for feature in record.features:
        for key in ['product', 'note', 'gene']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]
    return record.description

def has_protein(record):
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

        its1 = extract_region_status(description, "internal transcribed spacer 1")
        r5_8s = extract_region_status(description, "5.8S ribosomal RNA")
        its2 = extract_region_status(description, "internal transcribed spacer 2")
        lsu   = extract_region_status(description, "large subunit ribosomal RNA")

        function = extract_function(record, its1, r5_8s, its2, lsu)
        protein_flag = "Yes" if has_protein(record) else "-"

        # Save FASTA only if it doesn't already exist
        fasta_filename = os.path.join(output_folder, f"{accession}.fasta")
        save_fasta(record, fasta_filename)

        return accession, strain, its1, r5_8s, its2, lsu, function, protein_flag, description

    except Exception as e:
        return "Error", "-", "-", "-", "-", "-", "-", "-", str(e)

# Run the script
ids = fetch_accession_ids(search_term, max_results)
results = []

for idx, ncbi_id in enumerate(ids):
    row = fetch_metadata(ncbi_id)
    print(f"{idx+1}. {row[0]} | Protein: {row[7]} | Function: {row[6]}")
    results.append(row)
    time.sleep(0.4)

# Save to CSV
with open("pestalotiopsis_ITS_cleaned.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Accession", "Strain",
        "ITS1", "5.8S", "ITS2", "LSU",
        "Function", "Protein", "Full Description"
    ])
    writer.writerows(results)