import pandas as pd
import datetime
import matplotlib.pyplot as plt
import os
import numpy as np
import itertools
import re

timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# -------- CONFIGURATION --------
csv_file = r'C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Datasets\Results_accurate_Palladium.csv'
output_dir = 'Palladium-V2'
output_csv = os.path.join(output_dir, "Results_accurate_pv2.csv")
summary_file = f"summary-RU2{datetime.datetime.now():%Y%m%d_%H%M%S}.txt"
# --------------------------------

os.makedirs(output_dir, exist_ok=True)

# Detect if it's a Playground dataset
is_playground = "playground" in os.path.basename(csv_file).lower()

df = pd.read_csv(csv_file)
df.columns = df.columns.str.strip()
df.replace(["< LOD", "", " "], pd.NA, inplace=True)

if "Location" not in df.columns:
    raise ValueError("Missing 'Location' column in CSV!")

df["Location"] = df["Location"].replace({
    "Altadena - Hardware Store": "Altadena Hardware Store",
    "Temescal Gateway Park - Unburned": "Temescal Gateway Park",
    "Angeles National Forest - Main Trail": "Angeles National Forest",
    "Dodger's - Playground": "Dodger's Playground",
    "#3 - Tuna Canyon Park": "Tuna Canyon Park"
})

locations = df["Location"]

# Identify element columns
element_columns = [col for col in df.columns if not col.endswith("Err") and col not in ["File #","Location","DateTime","ElapsedTime","Cal Check"]]

for col in element_columns:
    df[col] = pd.to_numeric(df[col], errors='coerce') * 10000

for col in element_columns:
    df[f"{col}_zscore"] = (df[col] - df[col].mean()) / df[col].std()

df.to_csv(output_csv, index=False)
print(f"\n✅ Z-scores added and saved to: {output_csv}")

cols_lower_map = {col.lower(): col for col in df.columns}

element_colors = {
    "cu": "#d19563",
    "pb": "#7f7f7f",
    "pd": "#91d0e8",
    "y":  "#8dd3c7",
    "rb": "#bebada",
    "sr": "#fb8072",
    "ga": "#2ca02c"
}

element_full_names = {
    "Cu": "Copper",
    "Pd": "Palladium",
    "Y": "Yttrium",
    "Rb": "Rubidium",
    "Sr": "Strontium",
    "Pb": "Lead",
    "Zn": "Zinc",
    "Cr": "Chromium",
    "Sn": "Tin",
    "Sb": "Antimony",
    "Mo": "Molybdenum",
    "Ni": "Nickel",
    "Ga": "Gallium",
    "As": "Arsenic",
    "Ce": "Cerium"
}

fallback_colors = itertools.cycle(plt.cm.tab10.colors)
def get_color_for_element(element):
    el_lower = element.lower()
    if el_lower in element_colors:
        return element_colors[el_lower]
    else:
        new_color = next(fallback_colors)
        element_colors[el_lower] = new_color
        return new_color

def sanitize_filename(name):
    """Remove invalid Windows filename characters"""
    return re.sub(r'[<>:"/\\|?* ]', '_', name)

# ------------------ CHART FUNCTIONS ------------------

def plot_chart(elements, output_name_suffix):
    plt.figure(figsize=(10, 6))
    x = np.arange(len(locations))  # always use all locations for even spacing

    # Check if elements exist and have valid data
    for element in elements:
        error_col = element + " Err"
        if element not in df.columns or error_col not in df.columns:
            print(f"⚠️ Warning: Element '{element}' or its error column not found, skipping.")
            return
        values_check = pd.to_numeric(df[element], errors='coerce')
        if values_check.isna().all() or (values_check.fillna(0) == 0).all():
            print(f"⚠️ Warning: No valid data for elements {elements}, skipping chart.")
            return

    # Plot bars
    if len(elements) == 1:
        element = elements[0]
        values = pd.to_numeric(df[element], errors='coerce').fillna(0)
        errors_raw = pd.to_numeric(df[element + " Err"], errors='coerce') * 10000
        # Only keep error bars where the original value is not NaN
        errors = errors_raw.where(~df[element].isna(), 0)
        plt.bar(x, values, width=0.8, yerr=errors, capsize=3.5,
                color=get_color_for_element(element), edgecolor='black',
                error_kw={'elinewidth': 1})
    else:
        bottom = np.zeros(len(locations))
        for element in elements:
            values = pd.to_numeric(df[element], errors='coerce').fillna(0)
            plt.bar(x, values, bottom=bottom, width=0.8,
                    color=get_color_for_element(element), edgecolor='black')
            bottom += values

    plt.xlabel("Location")
    plt.ylabel("Concentration (ppm)")
    full_names = [f"{element_full_names.get(el, el)} ({el})" for el in elements]
    plt.title(", ".join(full_names) + " Concentrations in Los Angeles County Soils by Location", fontweight='bold')

    plt.xticks(x, locations, rotation=45, ha="right")

    ymax = plt.ylim()[1]
    plt.ylim(0, ymax * 1.2)

    # Draw headers and dotted lines
    if not is_playground:
        mid = len(locations) / 2
        plt.axvline(mid - 0.5, color='black', linestyle='--', linewidth=1.5)
        plt.text(mid / 2, ymax * 1.05, "Urban", ha='center', va='top')
        plt.text(3 * mid / 2, ymax * 1.05, "Recreational", ha='center', va='top')
    else:
        # Playground dataset: dotted lines and headers
        line1 = 2.5
        line2 = 7.5
        line3 = line2 + 5.0

        plt.axvline(line1, color='black', linestyle='--', linewidth=1.5)
        plt.axvline(line2, color='black', linestyle='--', linewidth=1.5)
        plt.axvline(line3, color='black', linestyle='--', linewidth=1.5)

        plt.text(line1 / 3, ymax * 1.15, "Controls", ha='center', va='top', fontsize=12)
        plt.text((line1 + line2) / 2, ymax * 1.15, "San Gabriel Valley\nPlaygrounds",
                 ha='center', va='top', fontsize=12)
        plt.text((line2 + line3) / 2, ymax * 1.15, "Burned\nBuildings",
                 ha='center', va='top', fontsize=12)
        plt.text((line3 + len(locations) - 0.1) / 2, ymax * 1.15,
                 "Wildfire Affected\nPark Areas", ha='center', va='top', fontsize=12)

    plt.tight_layout()

    # Save as PNG and SVG
    safe_name = sanitize_filename(output_name_suffix)
    png_path = os.path.join(output_dir, f"{safe_name}_ppm.png")
    plt.savefig(png_path, dpi=300)

    svg_path = os.path.join(output_dir, f"{safe_name}_ppm.svg")
    plt.savefig(svg_path, format="svg")

    plt.close()
    print(f"✅ Saved: {png_path} and {svg_path}")

# ------------------ RUN ------------------

summary_lines = []

print("\nAvailable elements in the dataset:")
for el in element_columns:
    print(f" - {el}")

# Ask whether to create single element charts
skip_singles_input = input("\nDo you want to generate single-element charts? (y/n): ").strip().lower()
skip_singles = (skip_singles_input == "n")

if not skip_singles:
    for element in element_columns:
        plot_chart([element], element.lower())

# Ask whether to create stacked charts
stacked_input = input("\nDo you want to generate stacked charts? (y/n): ").strip().lower()
if stacked_input == "y":
    print("\nEnter elements to stack (comma-separated, e.g. Cu,Pd,Pb):")
    user_input = input("Elements: ").strip()
    chosen_elements = [el.strip() for el in user_input.split(",") if el.strip()]

    # Validate chosen elements
    valid_elements = [el for el in chosen_elements if el in element_columns]
    invalid_elements = [el for el in chosen_elements if el not in element_columns]

    if invalid_elements:
        print(f"⚠️ Warning: These elements are not in dataset and will be skipped: {invalid_elements}")

    if valid_elements:
        plot_chart(valid_elements, "_stacked_" + "_".join(valid_elements).lower())
    else:
        print("❌ No valid elements selected. Skipping stacked chart.")

# Ask whether to create TOTAL SUM chart (all elements)
total_input = input("\nDo you want to generate a TOTAL SUM chart (all elements combined)? (y/n): ").strip().lower()
if total_input == "y":
    all_elements = []
    for base in element_columns:
        values = pd.to_numeric(df[base], errors='coerce')
        if not values.isna().all() and (values.fillna(0) != 0).any():
            all_elements.append(base)

    if all_elements:
        plot_chart(all_elements, "total_sum_all_elements")
    else:
        print("⚠️ Warning: No valid data found for any elements, skipping total sum chart.")

# Ask whether to create Zn + Pb total chart
zn_pb_input = input("\nDo you want to generate a Zn + Pb TOTAL chart? (y/n): ").strip().lower()
if zn_pb_input == "y":
    zn_pb = [el for el in ["Zn", "Pb"] if el in element_columns]
    if zn_pb:
        plot_chart(zn_pb, "total_sum_Zn_Pb")
    else:
        print("⚠️ Zn and/or Pb not found in dataset, skipping Zn+Pb total chart.")

# Save summary
summary_path = os.path.join(output_dir, summary_file)
with open(summary_path, 'w') as f:
    f.write("\n".join(summary_lines))

print(f"\n🎉 All charts saved in: {os.path.abspath(output_dir)}")
print(f"📄 Summary saved to: {os.path.abspath(summary_path)}")

