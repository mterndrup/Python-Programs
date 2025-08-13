from Bio import SeqIO

def merge_fasta_taxonomy():
    fasta_path = r"/BioTech/Bioinformatics/Tools/UNITE/QIIME/sh_qiime_release_19.02.2025/sh_refs_qiime_ver10_dynamic_19.02.2025.fasta"
    taxonomy_map_path = r"/BioTech/Bioinformatics/Tools/UNITE/QIIME/sh_qiime_release_19.02.2025/sh_taxonomy_qiime_ver10_dynamic_19.02.2025.txt"
    output_fasta_path = r"/BioTech/Bioinformatics/Tools/UNITE/QIIME/sh_qiime_release_19.02.2025/merged_sh_qiime_19.02.2025.fasta"

    # Read taxonomy mapping into dictionary {ID: taxonomy}
    tax_dict = {}
    with open(taxonomy_map_path, "r") as tax_file:
        for line in tax_file:
            if line.strip():
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    tax_dict[parts[0]] = parts[1]

    # Parse fasta and write merged fasta with taxonomy in header
    with open(output_fasta_path, "w") as out_fasta:
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq_id = record.id
            taxonomy = tax_dict.get(seq_id, "")
            if taxonomy:
                record.description = f"{seq_id} {taxonomy}"
            else:
                record.description = seq_id
            record.id = seq_id  # keep only seq ID, no spaces in header ID
            SeqIO.write(record, out_fasta, "fasta")

    print(f"Merged FASTA with taxonomy saved to:\n{output_fasta_path}")

if __name__ == "__main__":
    merge_fasta_taxonomy()
