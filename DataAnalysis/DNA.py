
from main_functions import list_files_in_current_folder

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
    
