import os
import re

# Root directory
base_dir = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Chaparral\Round2"

# Pattern for original .fastq files (excluding _trimmed, _filtered, etc.)
original_pattern = re.compile(r"^(.*?_pass_barcode\d{2}_[a-f0-9]+_[a-f0-9]+_[0-9]+)\.fastq$")

# Traverse all barcode folders
for folder in os.listdir(base_dir):
    barcode_path = os.path.join(base_dir, folder)
    raw_path = os.path.join(barcode_path, "Raw")

    if os.path.isdir(raw_path):
        output_lines = []
        combined_data = []
        output_filename = None

        print(f"\nProcessing folder: {folder}")

        for file in os.listdir(raw_path):
            match = original_pattern.match(file)
            if match:
                file_path = os.path.join(raw_path, file)
                with open(file_path, 'r') as f:
                    combined_data.append(f.read())
                output_lines.append(f"Included: {file}")
                print(f" - Included: {file}")
                if not output_filename:
                    prefix = match.group(1).rsplit("_", 1)[0]  # remove _0 or _1
                    output_filename = f"{prefix}_full.fastq"

        # Write combined FASTQ file if any originals were found
        if combined_data and output_filename:
            combined_path = os.path.join(barcode_path, output_filename)
            with open(combined_path, 'w') as out_f:
                out_f.write('\n'.join(combined_data))
            print(f" -> Combined FASTQ written to: {combined_path}")

            # Write log file
            log_path = os.path.join(barcode_path, f"{folder}_log.txt")
            with open(log_path, 'w') as log_f:
                log_f.write("\n".join(output_lines))
            print(f" -> Log written to: {log_path}")
        else:
            print(f" !! No original .fastq files found in: {raw_path}")
