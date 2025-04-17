import os, os.path, sys
import glob
from xml.etree import ElementTree

current_directory = os.getcwd()
xmls = os.path.join(current_directory, 'xmls')
folder_47 = os.path.join(xmls, '47')

def run():
    try:
        xml_files = glob.glob(folder_47 + "/*.xml")

        # Create a root element for the combined XML
        combined_root = ElementTree.Element("CombinedData")

        for xml_file in xml_files:
            data = ElementTree.parse(xml_file).getroot()
            # Append the data from each file to the combined root
            combined_root.append(data)

        # Create an ElementTree object with the combined root
        xml_element_tree = ElementTree.ElementTree(combined_root)

        # Write the combined XML to a file
        with open("results47.xml", "wb") as output_file:
            xml_element_tree.write(output_file)

        print("Combined XML written to results47.xml")
    except Exception as e:
        print(e)

# Call the function
run()