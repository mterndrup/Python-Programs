import os

current_directory = os.getcwd()
path = 0

def list_files_folders():
    try:
        path = int(input("Files (1) or Folders (2)?"))
        files = [f for f in os.listdir(current_directory) if os.path.isfile(f)]
        subfolders = [ f.path for f in os.scandir(current_directory) if f.is_dir() ]
        if path == 1:
            return files
        elif path == 2:
            return subfolders
    except Exception as e:
        print(f"An error occurred: {e}")
        return []

def navFolder():
    file_list = list_files_folders()
    dataframes = ''
    selected_filename = ''
    if not file_list:
        print("No files found in the current folder.")
    else:
        print("Files in the current folder:")
        for i, file in enumerate(file_list, start=1):
            print(f"{i}. {file}")

        selected_file_index = int(input("Enter the number for file/folder: "))

        if path == 1:
            try:
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
        elif path == 2:
            pass



