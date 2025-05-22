from Bio.Blast import NCBIWWW, NCBIXML
import re

folder = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\wildfirePlants-DNA-nanopore-sequence\fastq_pass\barcode34"
filename = "FBA73546_pass_barcode34_af591219_c9d44752_0_trimmed.fastq"
file_path = f"{folder}\\{filename}"


def assign_primer(hit_title, species):
    hit_title_lower = hit_title.lower()
    species_lower = species.lower()

    # Determine primer by hit description first
    if "16s" in hit_title_lower or "16s" in species_lower:
        return "16S"
    elif "18s" in hit_title_lower or "18s" in species_lower:
        return "18S"
    elif "5.8s" in hit_title_lower or "5.8s" in species_lower:
        # 5.8S is part of ITS regions, decide by species
        if any(x in species_lower for x in ["fungus", "fungal", "ascomycota", "basidiomycota", "yeast", "mushroom"]):
            return "FITS"
        elif any(x in species_lower for x in ["plant", "cactus", "angiosperm", "gymnosperm"]):
            return "PITS"
        else:
            return "ITS (unknown type)"
    else:
        # If no specific rRNA region mentioned, guess ITS type by species keywords
        if any(x in species_lower for x in ["fungus", "fungal", "ascomycota", "basidiomycota", "yeast", "mushroom"]):
            return "FITS"
        elif any(x in species_lower for x in ["plant", "cactus", "angiosperm", "gymnosperm"]):
            return "PITS"
        else:
            return "Unknown"


def parse_and_blast_fastq(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Skip malformed first read if only 3 lines (common in your files)
    if len(lines) >= 3 and lines[0].startswith('@') and not lines[2].startswith('+'):
        print("⚠️ Skipping malformed first read (only 3 lines).")
        lines = lines[3:]

    read_number = 1
    for i in range(0, len(lines), 4):
        try:
            header = lines[i].strip()
            sequence = lines[i + 1].strip()
            plus = lines[i + 2].strip()
            quality = lines[i + 3].strip()

            print(f"\n=== Read {read_number} ===")
            print(f"Sequence: {sequence}")

            # Run BLAST (nucleotide blast)
            print("Running BLAST...")
            result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            # Get the top hit info (best alignment)
            if blast_record.alignments:
                top_alignment = blast_record.alignments[0]
                top_hsp = top_alignment.hsps[0]
                hit_title = top_alignment.title

                # Extract species name inside brackets if possible
                species_match = re.search(r"\[([^\]]+)\]", hit_title)
                species = species_match.group(1) if species_match else "Unknown species"

                primer = assign_primer(hit_title, species)

                print(f"Top hit: {hit_title}")
                print(f"Species: {species}")
                print(f"Assigned primer: {primer}")
                print(f"Length: {top_alignment.length}")
                print(f"E-value: {top_hsp.expect}")
                print(f"Score: {top_hsp.score}")
                print(f"Query start: {top_hsp.query_start}, Subject start: {top_hsp.sbjct_start}")
                print(f"Identities: {top_hsp.identities}/{top_hsp.align_length}")
                print("------")
            else:
                print("No BLAST hits found.")

            read_number += 1

        except IndexError:
            print(f"⚠️ Incomplete read at index {i}. Skipping.")
            break
        except Exception as e:
            print(f"⚠️ Error blasting read {read_number}: {e}")
            read_number += 1
            continue


parse_and_blast_fastq(file_path)