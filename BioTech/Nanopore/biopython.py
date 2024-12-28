import os, Bio, re
from Bio.Blast import NCBIWWW
import chardet

email = ""
keys_folder = os.path.join(os.path.expanduser('~'), 'OneDrive', 'Documents', 'keys')
with open(os.path.join(keys_folder, 'ncbi.txt'), 'r') as file:
    content = file.read()
    email = re.search(r'email:\s*(\S+)', content).group(1)

current_directory = os.getcwd()
rawfiles_folder = os.path.join(current_directory, 'RawFiles')
file_path = os.path.join(rawfiles_folder, 'Split48.txt')
with open(file_path, 'r') as text_file:
    lines = text_file.readlines()

Bio.Blast.email = email

for i, line in enumerate(lines):
    try:
        if i > 1140:
            if line.strip():
                fasta_inc = ">con"+str(i)+"\n"+line
                print(fasta_inc)
                result_stream = NCBIWWW.qblast("blastn", "nt", fasta_inc)
                data = result_stream.read()

                with open(f"blast_results{i}.xml", "w") as xml_file:
                    xml_file.write(data)

                print("Printed xml file"+str(i)+"")
    except Exception as e:
        print(e)
