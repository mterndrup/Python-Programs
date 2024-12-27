import os, Bio, re
from Bio import Blast

email = ""
keys_folder = os.path.join(os.path.expanduser('~'), 'OneDrive', 'Documents', 'keys')
with open(os.path.join(keys_folder, 'ncbi.txt'), 'r') as file:
    content = file.read()
    email = re.search(r'email:\s*(\S+)', content).group(1)
    print(email)

current_directory = os.getcwd()
rawfiles_folder = os.path.join(current_directory, 'RawFiles')
print(current_directory)

Blast.tool
'biopython'
Blast.email = email