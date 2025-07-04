import os

input_file = '10CS22_orang-all-proteins.txt'

if not os.path.exists(input_file):
    print(f"Error: The input file '{input_file}' does not exist in the current directory.")
    print("Current directory contents:")
    print(os.listdir())
    exit(1)

with open(input_file, 'r', encoding='utf-8') as file:
    content = file.read()

lines = content.strip().splitlines()

def format_sequence(sequence, line_length):
    return "\n".join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))

with open("10CS22_orang-valid-proteins.txt", "w") as output_file:
    protein_sequence = []
    valid_protein = True
    header = ""
    started = False

    for line in lines:
        line = line.strip()  # Remove leading/trailing whitespace
        if ">" in line:
            if protein_sequence and valid_protein:
                formatted_sequence = format_sequence(''.join(protein_sequence), 60)
                print(f"{header}\n{formatted_sequence}\n")  # Ensure correct spacing in the console with 50 characters per line
                output_file.write(f"{header}\n{formatted_sequence}\n\n")  # Ensure correct spacing and line length in the file
            protein_sequence = []
            valid_protein = True
            header = line
            started = True
        elif started and valid_protein:
            if "*" in line:
                valid_protein = False
            else:
                protein_sequence.append(line)

    # Don't forget to write the last captured sequence if it's valid
    if protein_sequence and valid_protein:
        formatted_sequence = format_sequence(''.join(protein_sequence), 60)
        print(f"{header}\n{formatted_sequence}\n")  # Ensure correct spacing in the console with 50 characters per line
        output_file.write(f"{header}\n{formatted_sequence}\n\n")  # Ensure correct spacing and line length in the file
