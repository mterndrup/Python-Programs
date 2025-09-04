import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Hard-coded file path
FILE_PATH = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Datasets\Results_accurate_wildfire.csv"

# Mapping of element abbreviations to full names
element_full_names = {
    "As": "Arsenic",
    "Pb": "Lead",
    "Cd": "Cadmium",
    "Zn": "Zinc",
    "Cu": "Copper",
    "Ni": "Nickel",
    "Co": "Cobalt",
    "Cr": "Chromium",
    "Sb": "Antimony",
    "V": "Vanadium",
    "W": "Tungsten",
    "Mo": "Molybdenum",
    "Ga": "Gallium",
    "Te": "Tellurium",
    "Th": "Thorium",
    "Hf": "Hafnium",
    "Ce": "Cerium",
    "Ti": "Titanium",
    "Zr": "Zirconium",
    "Cl": "Chlorine",
    "Y": "Yttrium"
}

# Define risk categories
high_risk = ["As", "Co", "Cr", "Cu", "Ni", "Pb", "Sb", "V", "Zn"]
low_risk = ["Mo", "Ga", "Te", "Th", "Hf", "W", "Y"]

# Hard-coded top outlier sample name and extra info
TOP_OUTLIER_NAME = "Top Outlier: 822 E Mariposa St, Altadena"
TOP_OUTLIER_INFO = "Across from the burned hardware store"

def load_and_clean_csv(file_path):
    df = pd.read_csv(file_path)
    df = df.replace(r"< LOD", np.nan, regex=True)

    for col in df.columns:
        try:
            df[col] = df[col].astype(float)
        except ValueError:
            pass  # keep text columns like Location, GPS, etc.

    # Drop error columns and unnecessary metadata
    drop_cols = ["File #", "GPS", "ElapsedTime", "Cal Check"]
    df = df[[c for c in df.columns if not c.endswith("Err") and c not in drop_cols]]

    # Convert values to ppm
    for col in df.columns:
        df[col] = df[col].apply(lambda x: x * 10000 if pd.notnull(x) else np.nan)
    return df

def plot_boxplots(df, elements):
    selected_data = df[elements]

    # Determine colors based on risk
    colors = ['tab:blue' if el in high_risk else 'yellow' for el in elements]

    plt.figure(figsize=(max(12, len(elements)), 6))
    bp = plt.boxplot(
        [selected_data[el].dropna() for el in elements],
        patch_artist=True,
        showfliers=True
    )

    # Color boxes
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    # X-axis labels: abbreviation + full name
    new_labels = [f"{el}\n({element_full_names.get(el, '')})" for el in elements]
    plt.xticks(ticks=range(1, len(elements)+1), labels=new_labels, rotation=0, ha="center")

    # Vertical black line immediately to the right of Zn
    if "Zn" in elements:
        zn_index = elements.index("Zn") + 1
        plt.axvline(x=zn_index + 0.5, color='black', linestyle='--', linewidth=2)

    # Plot hard-coded top outlier
    for i, el in enumerate(elements):
        if el in df.columns:
            y_vals = df[el].dropna()
            if not y_vals.empty:
                y = y_vals.max()  # can be replaced with specific logic if known
                plt.plot(i + 1, y, 'o', color='red', markersize=8, label='_nolegend_')

    # Top-left legend: Top Outlier
    outlier_legend = plt.legend(
        handles=[plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10)],
        labels=[f"{TOP_OUTLIER_NAME}\n{TOP_OUTLIER_INFO}"],
        loc='upper left'
    )
    plt.gca().add_artist(outlier_legend)  # Keep this legend

    # Top-right legend: High Risk vs Low Risk
    risk_legend = plt.legend(
        handles=[
            plt.Line2D([0], [0], color='tab:blue', lw=10),
            plt.Line2D([0], [0], color='yellow', lw=10)
        ],
        labels=['High Risk', 'Low Risk'],
        loc='upper right'
    )
    plt.gca().add_artist(risk_legend)  # optional; sometimes plt.legend auto-adds

    plt.yscale("log")
    plt.ylabel("Concentration (ppm, log scale)")
    plt.title("Box Plot of Metals/Metalloids in Wildfire Affected Soils from LA County (n=32)")
    plt.grid(True, which="both", linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()


def main():
    df = load_and_clean_csv(FILE_PATH)

    print("\nAvailable elements:")
    print(", ".join(df.columns))

    elements = input("\nEnter comma-separated elements to plot (e.g., As,Co,Cr,Cu,Ni,Pb,Sb,V,Zn,Mo,Ga,Te,Th,Hf,W,Y): ").split(",")
    elements = [el.strip() for el in elements if el.strip() in df.columns]

    if not elements:
        print("No valid elements selected.")
        return

    plot_boxplots(df, elements)

if __name__ == "__main__":
    main()
