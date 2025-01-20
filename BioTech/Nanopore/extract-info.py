from collections import Counter, defaultdict

def extract_information(file_content):
    lines = file_content.split('\n')
    dna_sequence = ""
    genus_list = []
    species_list = []
    output = []
    current_kingdom = "Undetermined"

    for line in lines:
        if line.startswith("Line"):
            if dna_sequence:
                # Process the previous DNA sequence block
                most_common_kingdom = current_kingdom

                most_common_genus = Counter(genus_list).most_common(1)
                if most_common_genus and most_common_genus[0][1] > 1:
                    most_common_genus_name = f"{most_common_genus[0][0]} (Count: {most_common_genus[0][1]})"
                else:
                    most_common_genus_name = "Undetermined"

                output.append(f"{line_number}")
                output.append(f"{dna_sequence}")
                output.append(f"Most Common Genus: {most_common_genus_name}")
                output.append(f"Most Common Kingdom: {most_common_kingdom}")
                output.append("")

            # Reset for the next block
            line_number = line.split(": ", 1)[0]
            dna_sequence = line.split(": ", 1)[1]  # Extract the nucleotide sequence
            genus_list = []
        elif line.startswith("\t"):
            # Taxonomic data
            taxonomic_data = line.strip()
            if taxonomic_data.startswith("MAG:"):
                genus_list.append(taxonomic_data.split()[1])  # Take the next word after "MAG:"
            else:
                genus_list.append(taxonomic_data.split()[0])  # Assuming genus is the first word
            species_list.append(taxonomic_data)  # Full taxonomic data for species
        elif line.startswith("[{"):
            # Kingdom data in JSON format
            kingdom_data = eval(line)
            for item in kingdom_data:
                if isinstance(item, dict) and 'kingdom' in item:
                    current_kingdom = item['kingdom']


    # Count all species and sort by most common descending
    species_counter = Counter(species_list)
    sorted_species_counts = species_counter.most_common()

    # Organize species counts by genus
    genus_species_counts = defaultdict(list)
    for species, count in sorted_species_counts:
        genus_name = species.split()[0]
        genus_species_counts[genus_name].append((species, count))

    # Add species counts to the output organized by genus
    output.append("\nSpecies Counts by Genus:")
    for genus, species_counts in genus_species_counts.items():
        output.append(f"\nGenus: {genus}")
        for species, count in species_counts:
            output.append(f"  {species}: {count}")

    return "\n".join(output)

# Read the content of the original file
with open('output47-full.txt', 'r') as file:
    file_content = file.read()

# Extract information
extracted_info = extract_information(file_content)

# Write the result to a text file
with open('output47-extracted_info.txt', 'w') as output_file:
    output_file.write(extracted_info)

# Print the result to the terminal
print(extracted_info)

print("The extracted information has been saved to output47-extracted_info.txt")