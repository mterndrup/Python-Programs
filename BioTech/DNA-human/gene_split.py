import os

input_file = 'pGLO_plasmid_sequence.txt'

if not os.path.exists(input_file):
    print(f"Error: The input file '{input_file}' does not exist in the current directory.")
    print("Current directory contents:")
    print(os.listdir())
    exit(1)

with open(input_file, 'r', encoding='utf-8') as file:
    content = file.read()

lines = content.strip().split('ATG')
print(lines)