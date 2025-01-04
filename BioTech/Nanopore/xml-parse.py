import os, os.path, sys
import glob, re
from xml.etree import ElementTree
from Bio import Entrez
from collections import Counter

current_directory = os.getcwd()
xmls = os.path.join(current_directory, 'xmls')
folder_48 = os.path.join(xmls, '48')
rawfiles_folder = os.path.join(current_directory, 'RawFiles')

email = ""
keys_folder = os.path.join(os.path.expanduser('~'), 'OneDrive', 'Documents', 'keys')
with open(os.path.join(keys_folder, 'ncbi.txt'), 'r') as file:
    content = file.read()
    email = re.search(r'email:\s*(\S+)', content).group(1)

Entrez.email = email

class Logger:
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

sys.stdout = Logger("output48-full.txt")

def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid. we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace(' ', "+").strip()
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    return record['IdList'][0]

def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    return Entrez.read(search)

# Function to extract the numeric part of the file name
def extract_number(file_name):
    match = re.search(r'\d+', file_name)
    return int(match.group()) if match else float('inf')

kingdom_counter = []

def run():
    try:
        xml_files = glob.glob(folder_48 + "/*.xml")
        xml_files_sorted = sorted(xml_files, key=lambda x: extract_number(os.path.basename(x)))
        for xml_file in xml_files_sorted:
            file_name = os.path.basename(xml_file)
            number_match = re.search(r'\d+', file_name)
            if number_match:
                number = int(number_match.group())
                if number > -1:
                    print(f"Processing file number: {number}")

                    raw_file_path = os.path.join(rawfiles_folder, 'Split48.txt')
                    with open(raw_file_path, 'r') as raw_file:
                        lines = raw_file.readlines()
                        if number <= len(lines):
                            line_content = lines[number].strip()
                            print(f"Line {number} content: {line_content}")
                        else:
                            print(f"Line number {number} exceeds the total number of lines in the file.")

                    # Parse the XML file and extract <Hit_def> strings
                    tree = ElementTree.parse(xml_file)
                    root = tree.getroot()
                    hit_defs = [hit_def.text for hit_def in root.findall('.//Hit_def')]

                    # Extract and print the species (first two words or next two words if "PREDICTED:")
                    species_list = []
                    for hit_def in hit_defs:
                        words = hit_def.split()
                        if words[0] == "PREDICTED:":
                            species = ' '.join(words[1:3])
                        else:
                            species = ' '.join(words[:2])
                        species_list.append(species)

                    # Remove duplicates
                    unique_species_list = list(set(species_list))
                    taxid_list = []  # Initiate the lists to store the data to be parsed in
                    data_list = []
                    lineage_list = []
                    print('parsing taxonomic data...')  # message declaring the parser has begun

                    for species in unique_species_list:
                        try:
                            print('\t' + species)  # progress messages

                            taxid = get_tax_id(species)  # Apply your functions
                            data = get_tax_data(taxid)
                            #print(data)
                            lineage = {d['Rank']: d['ScientificName'] for d in data[0]['LineageEx'] if
                                       d['Rank'] in ['kingdom']}

                            taxid_list.append(taxid)  # Append the data to lists already initiated
                            data_list.append(data)
                            lineage_list.append(lineage)
                        except Exception as e:
                            print(e)
                    print(lineage_list)
                    try:
                        # Extract all kingdoms from the lineage_list
                        kingdoms = [list(lineage.values())[0] for lineage in lineage_list if lineage]

                        # Count the occurrences of each kingdom
                        kingdom_counts = Counter(kingdoms)

                        # Find the most common kingdom
                        most_common_kingdom = kingdom_counts.most_common(1)[0][0]

                        # Update kingdom_counter
                        kingdom_found = False
                        for kingdom in kingdom_counter:
                            if kingdom['name'] == most_common_kingdom:
                                kingdom['count'] += 1
                                kingdom_found = True
                                break
                        if not kingdom_found:
                            kingdom_counter.append({'name': most_common_kingdom, 'count': 1})

                        print(f'The most common kingdom is: {most_common_kingdom}')
                        print(f'Kingdom counter: {kingdom_counter}')
                    except Exception as e:
                        print(e)
                else:
                    print("Extracted number is not valid (must be greater than 0).")
            else:
                print("No number found in the file name.")
    except Exception as e:
        print(e)

run()

print(kingdom_counter)