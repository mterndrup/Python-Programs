import pandas as pd
import datetime
import matplotlib.pyplot as plt
import os
import numpy as np
import itertools

timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
# -------- CONFIGURATION --------
csv_file = r'C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Datasets\Results_accurate_Playgrounds-v2.csv'
output_dir = 'Palladium-PlayV2'
output_csv = os.path.join(output_dir, "Results_accurate_playv2.csv")
summary_file = f"summary-RU2{datetime.datetime.now():%Y%m%d_%H%M%S}.txt"
# --------------------------------

os.makedirs(output_dir, exist_ok=True)

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
    "pd": "#91d0e8",
    "y":  "#8dd3c7",
    "rb": "#bebada",
    "sr": "#fb8072"
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

# ------------------ CHART FUNCTIONS ------------------

def plot_chart(elements, output_name_suffix):
    plt.figure(figsize=(10, 6))
    x = np.arange(len(locations))

    for element in elements:
        error_col = element + " Err"
        if element not in df.columns or error_col not in df.columns:
            print(f"⚠️ Warning: Element '{element}' or its error column not found, skipping.")
            return
        values = pd.to_numeric(df[element], errors='coerce') * 10000
        if values.isna().all() or (values.fillna(0) == 0).all():
            print(f"⚠️ Warning: No valid data for elements {elements}, skipping chart.")
            return

    if len(elements) == 1:
        element = elements[0]
        values = pd.to_numeric(df[element], errors='coerce')
        errors = pd.to_numeric(df[element + " Err"], errors='coerce') * 10000
        plt.bar(x, values, width=0.8, yerr=errors, capsize=3.5,
                color=get_color_for_element(element), edgecolor='black', error_kw={'elinewidth': 1})
    else:
        bottom = np.zeros(len(locations))
        for element in elements:
            values = pd.to_numeric(df[element], errors='coerce') * 10000
            plt.bar(x, values, bottom=bottom, width=0.8,
                    color=get_color_for_element(element), edgecolor='black')
            bottom += np.nan_to_num(values)

    plt.xlabel("Location")
    plt.ylabel("Concentration (ppm)")
    plt.title(", ".join(elements) + " Concentrations by Location", fontweight='bold')
    plt.xticks(x, locations, rotation=45, ha="right")
    plt.axvline(len(locations)/2-0.5, color='black', linestyle='--', linewidth=1.5)
    plt.text(len(locations)/4, plt.ylim()[1]*0.95, "Urban", ha='center', va='top')
    plt.text(3*len(locations)/4, plt.ylim()[1]*0.95, "Recreational", ha='center', va='top')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{output_name_suffix}_ppm.png"))
    plt.close()
    print(f"✅ Saved: {output_name_suffix}_ppm.png")

def plot_combined_avg_chart(elements):
    if len(elements) != 2:
        print("⚠️ Please provide exactly two elements.")
        return

    mid = len(df)//2
    values = []
    errors = []

    for el in elements:
        vals = pd.to_numeric(df[el], errors='coerce')
        errs = pd.to_numeric(df[el + " Err"], errors='coerce') * 10000
        urban_val = np.nanmean(vals[:mid])
        rec_val = np.nanmean(vals[mid:])
        urban_err = np.nanmean(errs[:mid])
        rec_err = np.nanmean(errs[mid:])
        values.append([urban_val, rec_val])
        errors.append([urban_err, rec_err])

    x = np.arange(2)
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 6))

    legend_labels = [f"{element_full_names.get(el, el)} ({el})" for el in elements]

    ax.bar(x - width/2, [values[0][0], values[0][1]], width, yerr=[errors[0][0], errors[0][1]],
           capsize=5, label=legend_labels[0], color=get_color_for_element(elements[0]), edgecolor='black')
    ax.bar(x + width/2, [values[1][0], values[1][1]], width, yerr=[errors[1][0], errors[1][1]],
           capsize=5, label=legend_labels[1], color=get_color_for_element(elements[1]), edgecolor='black')

    ax.set_xticks(x)
    ax.set_xticklabels(['Urban', 'Recreational'])
    ax.set_ylabel("Average Concentration (ppm)")
    ax.set_title(f"Urban vs Recreational Averages: {legend_labels[0]} and {legend_labels[1]}", fontweight='bold')
    ax.legend()
    plt.tight_layout()

    chart_file = os.path.join(output_dir, f"{elements[0]}_{elements[1]}_urban_recreational_avg.png")
    plt.savefig(chart_file)
    plt.close()
    print(f"✅ Saved combined average chart with error bars: {chart_file}")

def plot_total_concentration_stacked(elements):
    if len(elements) != 2:
        print("⚠️ Please provide exactly two elements.")
        return

    fig_width = max(12, len(df) * 0.3)
    fig, ax = plt.subplots(figsize=(fig_width, 8))
    x = np.arange(len(elements))
    width = 0.6

    mid = len(df) // 2
    urban_samples = locations[:mid].tolist()
    rec_samples = locations[mid:].tolist()

    bottom = np.zeros(len(elements))

    # Use tab20 discrete colors
    cmap = plt.colormaps['tab20'].colors
    urban_colors = [cmap[i % 10] for i in range(mid)]
    rec_colors = [cmap[10 + i % 10] for i in range(len(df)-mid)]

    # Track cumulative heights to add horizontal line
    cumulative_urban = np.zeros(len(elements))

    for i, sample_idx in enumerate(range(len(df))):
        sample_values = []
        for j, el in enumerate(elements):
            val = pd.to_numeric(df.at[sample_idx, el], errors='coerce')
            sample_values.append(val if not np.isnan(val) else 0)

        color = urban_colors[i] if i < mid else rec_colors[i - mid]
        ax.bar(x, sample_values, width, bottom=bottom, color=color, edgecolor='black')
        bottom += sample_values

        if i == mid - 1:
            cumulative_urban = bottom.copy()  # cumulative height of last urban sample

    # Draw horizontal line slightly past each bar
    for j in range(len(elements)):
        line_start = x[j] - width * 0.7
        line_end = x[j] + width * 0.7
        ax.hlines(y=cumulative_urban[j], xmin=line_start, xmax=line_end,
                  colors='black', linestyles='--', linewidth=1.5)

    ax.set_xticks(x)
    element_labels = [f"{element_full_names.get(el, el)} ({el})" for el in elements]
    ax.set_xticklabels(element_labels)
    ax.set_ylabel("Concentration (ppm)")
    ax.set_title(f"Stacked Concentrations Across Samples: {element_labels[0]} and {element_labels[1]}", fontweight='bold')

    # Legends
    handles_urban = [plt.Rectangle((0,0),1,1,color=urban_colors[i]) for i in range(mid)]
    handles_rec = [plt.Rectangle((0,0),1,1,color=rec_colors[i]) for i in range(len(df)-mid)]

    urban_legend = ax.legend(handles_urban, urban_samples, title="Urban Samples",
                             bbox_to_anchor=(1.05, 1), loc='upper left')
    rec_legend = ax.legend(handles_rec, rec_samples, title="Recreational Samples",
                            bbox_to_anchor=(1.05, 0.5), loc='upper left')
    ax.add_artist(urban_legend)

    plt.tight_layout()
    chart_file = os.path.join(output_dir, f"{elements[0]}_{elements[1]}_stacked_samples.png")
    plt.savefig(chart_file, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved stacked total concentration chart: {chart_file}")


# ------------------ RUN ------------------

summary_lines = []

print("\nAvailable elements in the dataset:")
for el in element_columns:
    print(f" - {el}")

# Prompt for Urban/Recreational chart
combo_avg_input = input("\nDo you want a combined Urban/Recreational chart for two elements? Enter two elements separated by '+', or 'no' to skip: ").strip().lower()
if combo_avg_input != 'no' and combo_avg_input:
    elements_input = [el.strip() for el in combo_avg_input.split("+")]
    elements = []
    for el in elements_input:
        if el in cols_lower_map:
            elements.append(cols_lower_map[el])
        else:
            print(f"⚠️ Warning: Element '{el}' not found in data columns.")
    if len(elements) == 2:
        plot_combined_avg_chart(elements)
    else:
        print("⚠️ Error: Please provide exactly two valid elements.")

# Prompt for stacked total concentration chart
total_input = input("\nDo you want a stacked total concentration chart for two elements? Enter two elements separated by '+', or 'no' to skip: ").strip().lower()
if total_input != 'no' and total_input:
    elements_input = [el.strip() for el in total_input.split("+")]
    elements = []
    for el in elements_input:
        if el in cols_lower_map:
            elements.append(cols_lower_map[el])
        else:
            print(f"⚠️ Warning: Element '{el}' not found in data columns.")
    if len(elements) == 2:
        plot_total_concentration_stacked(elements)
    else:
        print("⚠️ Error: Please provide exactly two valid elements.")

# Plot single element charts
skip_singles = False  # you can still add a prompt to skip if needed
if not skip_singles:
    for element in element_columns:
        plot_chart([element], element.lower())

# Save summary
summary_path = os.path.join(output_dir, summary_file)
with open(summary_path, 'w') as f:
    f.write("\n".join(summary_lines))

print(f"\n🎉 All charts saved in: {os.path.abspath(output_dir)}")
print(f"📄 Summary saved to: {os.path.abspath(summary_path)}")
