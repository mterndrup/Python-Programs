import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import itertools

# -------- CONFIGURATION --------
csv_file = r'C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Results_accurate_RU2.csv'
output_dir = 'RU2'
summary_file = 'summary-RU2.txt'
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

locations = df["Location"]

# Identify element columns and their corresponding error columns
element_columns = []
for col in df.columns:
    if col.endswith(" Err"):
        base_col = col.replace(" Err", "")
        if base_col in df.columns:
            element_columns.append((base_col, col))

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
    for base, _ in element_columns:
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
        values = pd.to_numeric(df[element], errors='coerce') * 10000
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
    plt.title(f"Elements: {', '.join(elements)} Concentration by Location")
    plt.xticks(x, locations, rotation=45, ha="right")
    plt.tight_layout()

    chart_file = os.path.join(output_dir, f"{output_name_suffix}_ppm_with_error.png")
    plt.savefig(chart_file)
    plt.close()
    print(f"✅ Saved: {chart_file}")

# Plot individual element charts if not skipping
if not skip_singles:
    for element, error_col in element_columns:
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
