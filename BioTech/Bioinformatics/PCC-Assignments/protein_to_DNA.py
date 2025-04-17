import os

protein_file = '10CS22_orang-valid-proteins.txt'
DNA_file = '10CS22_orang.fasta'
output_file = '10CS22_orang-valid-DNA_sequences.txt'

def read_fasta(file_path):
    if not os.path.exists(file_path):
        print(f"Error: The input file '{file_path}' does not exist in the current directory.")
        print("Current directory contents:")
        print(os.listdir())
        exit(1)

    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()

    lines = content.strip().splitlines()
    data = {}
    current_identifier = None

    for line in lines:
        if line.startswith('>'):
            current_identifier = line[1:].strip()  # Remove '>' and any extra spaces
            data[current_identifier] = ''
        elif current_identifier:  # Ensure current_identifier is set
            data[current_identifier] += line.strip()

    return data

def split_sequence(sequence, length=50):
    return [sequence[i:i+length] for i in range(0, len(sequence), length)]

# Read protein and DNA files
proteins = read_fasta(protein_file)
dna_sequences = read_fasta(DNA_file)

# Write the matching DNA sequences to a new file with lines of max 60 characters
count = 0
with open(output_file, 'w', encoding='utf-8') as file:
    for identifier in proteins.keys():
        # Extract the part after the "rf 1 " in protein identifier to match DNA identifier
        dna_identifier = identifier.split()[2] if len(identifier.split()) > 2 else identifier
        if dna_identifier in dna_sequences:
            file.write(f'>{dna_identifier}\n')
            split_seq = split_sequence(dna_sequences[dna_identifier])
            file.write('\n'.join(split_seq) + '\n\n')  # Add an extra newline for spacing
            count += 1

print(f"Filtered DNA sequences have been written to '{output_file}'.")
print(f"Total number of sequences: {count}")