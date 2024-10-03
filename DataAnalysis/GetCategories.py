import os
import pandas as pd
import numpy as np
import csv

from main_functions import list_files_in_current_folder

# Set display options to avoid truncation
pd.set_option("display.max_colwidth", None)  # Display full column content
pd.set_option("display.max_rows", None)  # Display all rows

script_dir = os.path.dirname(os.path.realpath(__file__))
df_dict = ''


def read_all_sheets_from_file(filename) -> dict:
    try:
        path = os.path.join(script_dir, filename)
        _, ext = os.path.splitext(filename.lower())
        
        if ext == '.xlsx' or ext == '.xls':
            xls = pd.ExcelFile(path)
            df_dict = {sheet_name: xls.parse(sheet_name) for sheet_name in xls.sheet_names}
            df = pd.read_excel(filename, sheet_name=0)
        elif ext == '.csv':
            df = pd.read_csv(path)
            df_dict = {'Sheet1': df}  # You can customize the sheet name as needed
        else:
            print(f"Unsupported file format: {ext}. Only Excel (.xlsx, .xls) and CSV (.csv) files are supported.")
            return {}
        
        return df_dict
    except FileNotFoundError:
        print(f"File '{filename}' not found. Make sure it's in the same folder as the Python script.")
        return {}

def get_unique_values_from_csv(file_path, column1_name, column2_name):
    unique_values = set()
    try:
        dfi = pd.read_csv(file_path)
        dfi[column1_name] = pd.to_numeric(dfi[column1_name], errors='coerce')  # Convert column1 to numeric
        df_sorted = dfi.sort_values(by=column1_name, ignore_index=True)  # Sort by column1
        unique_values = set(zip(df_sorted[column1_name], df_sorted[column2_name]))
        return list(unique_values)
    except Exception as e:
        print(f"Error reading the file: {e}")
        return None
    

if __name__ == "__main__":
    file_list = list_files_in_current_folder()
    dataframes = ''
    selected_filename = ''
    if not file_list:
        print("No files found in the current folder.")
    else:
        print("Files in the current folder:")
        for i, file in enumerate(file_list, start=1):
            print(f"{i}. {file}")
        selected_file_index = input("Enter the number corresponding to the Excel file you want to read: ")
        try:
            selected_file_index = int(selected_file_index)
            if 1 <= selected_file_index <= len(file_list):
                selected_filename = file_list[selected_file_index - 1]
                dataframes = read_all_sheets_from_file(selected_filename)
                if not dataframes:
                    print(f"No data loaded from '{selected_filename}'.")
                else:
                    print(f"Data loaded successfully from '{selected_filename}'!")
            else:
                print("Invalid selection. Please choose a valid file number.")
        except ValueError:
            print("Invalid input. Please enter a valid file number.")
    if not dataframes:
        print("No data loaded due to missing Excel file.")
    else:
        #print("Data loaded successfully!")
        #print("What information do you want to Analyze?")
        #print("1. Category")
        #print("2. Subcategory")
        #sort_index = input("Type the number here: ")
        #print(sort_index)

        # Start with aggregating the data by city
        #df = dataframes['PART_I_AND_II_CRIMES_YTD_0']
        #df.drop_duplicates(subset=['Stat Code'], inplace=True)
        #all_categories = df[["Stat Code","Stat Code Desc"]].to_numpy()

        column_to_extract1 = 'Stat Code'
        column_to_extract2 = 'Stat Code Desc'
        unique_values_list = get_unique_values_from_csv(selected_filename, column_to_extract1, column_to_extract2)
        sorted_unique_values = sorted(unique_values_list, key=lambda x: float(x[0]))

        if sorted_unique_values:
            print("Unique values from both columns (sorted by value1):")
            for value1, value2 in sorted_unique_values:
                print(f"{value1} - {value2}")
        else:
            print("Error occurred while processing the file.")

        

        # Sort the unique category values alphabetically
        #sorted_categories = sorted(unique_categories)
        #sorted_categories = all_categories[np.argsort(all_categories[:, 0])]

        #print("Unique category values:")
        #for category in sorted_categories:
            #print(category)

