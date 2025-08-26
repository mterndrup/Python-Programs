import pandas as pd
import datetime
import matplotlib.pyplot as plt
import os
import numpy as np
import itertools

timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
# -------- CONFIGURATION --------
csv_file = r'C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Results_accurate_Palladium.csv'
output_dir = 'Palladium-Final'
output_csv = os.path.join(output_dir, "Results_accurate_Palladium_v2.csv")
summary_file = f"summary-RU2{datetime.datetime.now():%Y%m%d_%H%M%S}.txt"
# --------------------------------

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Read the CSV (comma-separated)
df = pd.read_csv(csv_file)

# Clean column names
df.columns = df.columns.str.strip()

# Replace '< LOD' and blanks with NaN
df.replace(["< LOD", "", " "], pd.NA, inplace=True)

# Check for 'Location'
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


# Identify element columns and their corresponding error columns
element_columns = [col for col in df.columns if not col.endswith("Err") and col not in ["File #","Location","DateTime","ElapsedTime","Cal Check"]]

# Convert all element columns to numeric (coerce errors)
for col in element_columns:
    df[col] = pd.to_numeric(df[col], errors='coerce') * 10000  # convert to ppm

# Calculate z-scores for each element
for col in element_columns:
    mean_val = df[col].mean()
    std_val = df[col].std()
    z_col_name = f"{col}_zscore"
    df[z_col_name] = (df[col] - mean_val) / std_val

# Save new CSV
df.to_csv(output_csv, index=False)
print(f"\n✅ Z-scores added and saved to: {output_csv}")

# Create mapping lowercase -> actual DataFrame column name
cols_lower_map = {col.lower(): col for col in df.columns}

# Map for bar colors (keys are lowercase) — add your preferred colors here
element_colors = {
    "cu": "#d19563",   # copper
    "pd": "#91d0e8",   # palladium
    "y":  "#8dd3c7",   # yttrium
    "rb": "#bebada",   # rubidium
    "sr": "#fb8072"    # strontium
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
    "Ce": "Cerium",
    # add others if needed...
}

# Fallback color cycle for unlisted elements
fallback_colors = itertools.cycle(plt.cm.tab10.colors)

def get_color_for_element(element):
    """Return a unique color for each element."""
    el_lower = element.lower()
    if el_lower in element_colors:
        return element_colors[el_lower]
    else:
        new_color = next(fallback_colors)
        element_colors[el_lower] = new_color  # store for consistency
        return new_color

# Ask user about combined charts
combo_input = input(
    "Enter elements to combine in charts (e.g. 'cu+pd, pd+pb+cu'), "
    "'total' for all elements, or 'skip' to skip singles: "
).strip()

combo_groups = []
skip_singles = False
add_total = False

if combo_input:
    parts = [p.strip().lower() for p in combo_input.split(",") if p.strip()]
    if "skip" in parts:
        skip_singles = True
        parts.remove("skip")
    if "total" in parts:
        add_total = True
        parts.remove("total")

    for combo in parts:
        combo_elements = []
        for el in combo.split("+"):
            el_clean = el.strip().lower()
            if el_clean in cols_lower_map:
                combo_elements.append(cols_lower_map[el_clean])
            else:
                print(f"⚠️ Warning: Element '{el_clean}' not found in data columns.")
        if combo_elements:
            combo_groups.append(combo_elements)

# If user requested 'total', build it after filtering out empty elements
if add_total:
    all_elements = []
    for base in element_columns:   # not for base, _
        values = pd.to_numeric(df[base], errors='coerce') * 10000
        if not values.isna().all() and (values.fillna(0) != 0).any():
            all_elements.append(base)
    if all_elements:
        combo_groups.append(all_elements)
    else:
        print("⚠️ Warning: No valid data found for any elements.")

summary_lines = []

def plot_chart(elements, output_name_suffix):
    """Plot chart for given elements list."""
    plt.figure(figsize=(10, 6))
    x = np.arange(len(locations))

    # Check if any element has all NaN or zero values; skip chart if so
    for element in elements:
        error_col = element + " Err"
        if element not in df.columns or error_col not in df.columns:
            print(f"⚠️ Warning: Element '{element}' or its error column not found, skipping.")
            return
        values = pd.to_numeric(df[element], errors='coerce') * 10000
        if values.isna().all() or (values.fillna(0) == 0).all():
            print(f"⚠️ Warning: No valid data to plot for elements {elements}, skipping chart.")
            return

    if len(elements) == 1:
        # Single element with error bars
        element = elements[0]
        error_col = element + " Err"
        values = pd.to_numeric(df[element], errors='coerce')
        errors = pd.to_numeric(df[error_col], errors='coerce') * 10000
        color = get_color_for_element(element)
        plt.bar(x, values, width=0.8, yerr=errors, capsize=3.5,
                color=color, edgecolor='black', error_kw={'elinewidth': 1})

        avg_ppm = np.nanmean(values)
        avg_error = np.nanmean(errors)
        efficiency_error_rate = (avg_error / avg_ppm) * 100 if avg_ppm and not np.isnan(avg_ppm) else np.nan
        summary_lines.append(
            f"{element}: Average ppm = {avg_ppm:.3f}, Average error = {avg_error:.3f}, "
            f"Efficiency error rate = {efficiency_error_rate:.2f}%"
        )

    else:
        # Stacked bar chart for combos (no error bars)
        bottom = np.zeros(len(locations))
        for element in elements:
            values = pd.to_numeric(df[element], errors='coerce') * 10000
            color = get_color_for_element(element)
            plt.bar(x, values, bottom=bottom, width=0.8,
                    color=color, edgecolor='black')
            bottom += np.nan_to_num(values)

            avg_ppm = np.nanmean(values)
            summary_lines.append(
                f"{element}: Average ppm = {avg_ppm:.3f} (error bars omitted in stacked combo charts)"
            )

    plt.xlabel("Location")
    plt.ylabel("Concentration (ppm)")
    # Convert symbols to "Name (Symbol)" form
    pretty_elements = []
    for el in elements:
        symbol = el.strip()
        name = element_full_names.get(symbol, symbol)  # fallback to symbol if missing
        pretty_elements.append(f"{name} ({symbol})")

    plt.title(f"{', '.join(pretty_elements)} Concentrations in Los Angeles County Soils (2025) by Location",fontweight='bold')
    plt.xticks(x, locations, rotation=45, ha="right")
    plt.tight_layout()

    # Draw vertical line in the middle
    middle_x = len(locations) / 2 - 0.5  # subtract 0.5 to align with bar edges
    plt.axvline(x=middle_x, color='black', linestyle='--', linewidth=1.5)

    # Add labels for each half
    plt.text(middle_x / 2, plt.ylim()[1] * 0.95, "Urban", ha='center', va='top', fontsize=12)
    plt.text(middle_x + middle_x / 2 + 0.5, plt.ylim()[1] * 0.95, "Recreational", ha='center', va='top', fontsize=12)

    chart_file = os.path.join(output_dir, f"{output_name_suffix}_ppm_with_error.png")
    plt.savefig(chart_file)
    plt.close()
    print(f"✅ Saved: {chart_file}")

# Plot individual element charts if not skipping
if not skip_singles:
    for element in element_columns:  # not for element, error_col
        plot_chart([element], element.lower())

# Plot combined charts if any
for combo in combo_groups:
    combo_name = "_".join([e.lower() for e in combo])
    plot_chart(combo, combo_name)

# Save summary in the output directory
summary_path = os.path.join(output_dir, summary_file)
with open(summary_path, 'w') as f:
    f.write("\n".join(summary_lines))

print(f"\n🎉 All charts saved in: {os.path.abspath(output_dir)}")
print(f"📄 Summary saved to: {os.path.abspath(summary_path)}")
