import csv
import pandas as pd
from dateutil import parser

# --- Ask user for date ---
target_date_str = input("Enter a date to process samples from (any common format): ").strip()

try:
    target_date = parser.parse(target_date_str, dayfirst=False).date()
except ValueError:
    print("Invalid date format. Please enter a recognizable date.")
    exit()

# --- File paths ---
file_path = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Results.csv"
excel_output = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Results_output.xlsx"

# --- Helpers ---
def is_header_row(row):
    return any(h in row for h in ['File #', 'DateTime', 'Operator'])

def is_data_row(row):
    for cell in row:
        cell = cell.strip()
        if cell == "<LOD":
            return True
        try:
            float(cell)
            return True
        except:
            continue
    return False

# --- Parse file ---
samples = []
with open(file_path, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rows = list(reader)

i = 0
while i < len(rows):
    row = rows[i]

    if not any(cell.strip() for cell in row):
        i += 1
        continue

    if is_header_row(row):
        header_row = [h.strip() for h in row]

        if i + 1 < len(rows) and is_data_row(rows[i + 1]):
            value_row = rows[i + 1]
            sample = {}
            for header, value in zip(header_row, value_row):
                value = value.strip()
                if header:
                    if value == "<LOD":
                        sample[header] = None
                    else:
                        try:
                            sample[header] = float(value)
                        except ValueError:
                            sample[header] = value
            # Filter by date
            dt_val = sample.get("DateTime", "")
            sample_date_str = str(dt_val).split()[0] if dt_val else ""
            try:
                sample_date = parser.parse(sample_date_str).date()
                if sample_date == target_date:
                    samples.append(sample)
            except:
                pass
            i += 2
        else:
            i += 1
    else:
        i += 1

if not samples:
    print(f"No samples found for date: {target_date}")
    exit()

# --- Detect metadata columns dynamically ---
first_err_idx = next((idx for idx, h in enumerate(header_row) if h.endswith("Err")), None)
if first_err_idx is None:
    print("No 'Err' columns found, cannot determine elements vs metadata.")
    exit()

first_element_idx = first_err_idx - 1  # first element column
metadata_cols = header_row[:first_element_idx]  # everything before first element is metadata

# --- Process each sample ---
sorted_samples = []
all_elements = set()

for sample in samples:
    metadata_keys = [k for k in sample.keys() if k in metadata_cols]
    element_error_keys = [k for k in sample.keys() if k not in metadata_cols]

    # Process elements/errors in header order
    element_error_keys_sorted = [k for k in header_row if k in element_error_keys]

    elements = {}
    skip_next = False
    for i, key in enumerate(element_error_keys_sorted):
        if skip_next:
            skip_next = False
            continue
        if key.endswith("Err"):
            continue
        val = sample.get(key)
        err = None
        if i + 1 < len(element_error_keys_sorted) and element_error_keys_sorted[i + 1] == f"{key} Err":
            err = sample.get(element_error_keys_sorted[i + 1])
            skip_next = True
        elements[key] = {"Val": val, "Err": err}
        all_elements.add(key)

    # Sort elements alphabetically
    sorted_element_names = sorted(elements.keys(), key=lambda x: x.lower())

    flat_sample = {}
    for mk in metadata_cols:
        flat_sample[mk] = sample.get(mk)
    for elem_name in sorted_element_names:
        flat_sample[elem_name] = elements[elem_name]["Val"]
        if elements[elem_name]["Err"] is not None:
            flat_sample[f"{elem_name} Err"] = elements[elem_name]["Err"]

    sorted_samples.append(flat_sample)

# --- Final columns ---
final_metadata_cols = [col for col in metadata_cols if any(col in s for s in sorted_samples)]
sorted_all_elements = sorted(all_elements, key=lambda x: x.lower())

final_columns = final_metadata_cols.copy()
for elem_name in sorted_all_elements:
    final_columns.append(elem_name)
    err_col = f"{elem_name} Err"
    if any(err_col in s for s in sorted_samples):
        final_columns.append(err_col)

# --- Export to Excel ---
df = pd.DataFrame(sorted_samples, columns=final_columns)
df.to_excel(excel_output, index=False)
print(f"\nData exported to Excel:\n{excel_output}")

