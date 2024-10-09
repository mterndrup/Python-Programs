import os
from main_functions import navFolder
input_file = 'ManuSporny-genome.txt'

if not os.path.exists(input_file):
    print(f"Error: The input file '{input_file}' does not exist in the current directory.")
    print("Current directory contents:")
    print(os.listdir())
    exit(1)

with open(input_file, 'r', encoding='utf-8') as file:
    content = file.read()

lines = content.strip().split('\n')
chromosomes = [[] for i in range(25)]


for line in lines:
    parts = line.split()
    if len(parts) == 4:
        rsid, chrm, pos, geno = parts
        chrmSection = [parts]
        if chrm == 'X':
            chromosomes[int(23)].append(chrmSection)
        elif chrm == 'Y':
            chromosomes[int(24)].append(chrmSection)
        elif chrm == 'MT':
            chromosomes[0].append(chrmSection)
        elif int(chrm) <=22:
            chromosomes[int(chrm)].append(chrmSection)
        print(".", end='')
            
            

print(len(chromosomes))
for i, item in enumerate(chromosomes):
    with open("Output"+str(i)+".txt", "w") as text_file:
        for z, ch in enumerate(item):
            output = str(ch)
            output = output.replace('[','')
            output = output.replace(']','')
            output = output.replace('\'','')
            output = output.replace('\'','')
            text_file.write(output+"\n")
        print(len(chromosomes[i]))

    

