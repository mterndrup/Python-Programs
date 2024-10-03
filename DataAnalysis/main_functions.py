import os

current_directory = os.getcwd()

def list_files_in_current_folder():
    try:
        path = int(input("Files (1) or Folders (2)?"))
        files = [f for f in os.listdir(current_directory) if os.path.isfile(f)]
        subfolders = [ f.path for f in os.scandir(current_directory) if f.is_dir() ]
        return files
    except Exception as e:
        print(f"An error occurred: {e}")
        return []




