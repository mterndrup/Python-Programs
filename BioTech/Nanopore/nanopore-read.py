import os

input_file = 'FAY23443_pass_barcode47_b0113415_c07d72eb_0.fastq'

if not os.path.exists(input_file):
    print(f"Error: The input file '{input_file}' does not exist in the current directory.")
    print("Current directory contents:")
    print(os.listdir())
    exit(1)

with open(input_file, 'r', encoding='utf-8') as file:
    content = file.read()

lines = content.strip().splitlines()
with open("basepairmedium47.txt", "w") as text_file:
    for i, line in enumerate(lines):
        try:
            output = str(lines[i-1])
            print(lines[i-1])
            text_file.write(f">con{i} \n {output}\n\n")
        except Exception as e:
            print(e)
