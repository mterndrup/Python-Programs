import json
import csv
import os
from datetime import datetime
import matplotlib.pyplot as plt
from pdz_tool import PDZTool

# --- File paths ---
pdz_file = r"C:\Users\ketgl\OneDrive\Documents\GitHub\HTML-Web\XRF-Processor\00357-GeoExploration.pdz"
json_out = r"C:\Users\ketgl\OneDrive\Documents\GitHub\HTML-Web\XRF-Processor\output.json"
csv_out = r"C:\Users\ketgl\OneDrive\Documents\GitHub\HTML-Web\XRF-Processor\output.csv"

# Timestamped CSV to avoid lock issues
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
csv_out = csv_out.replace(".csv", f"_{timestamp}.csv")

# --- Initialize PDZ tool ---
pdz_tool = PDZTool(pdz_file)

# --- Parse file ---
parsed_data = pdz_tool.parse()

# --- Save JSON ---
with open(json_out, "w", encoding="utf-8") as f:
    json.dump(parsed_data, f, indent=2)
print(f"✅ Saved JSON to {json_out}")

# --- Save CSV (flattened) ---
try:
    with open(csv_out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Record Name", "Field", "Value"])
        for record_name, record_data in parsed_data.items():
            if isinstance(record_data, dict):
                for k, v in record_data.items():
                    writer.writerow([record_name, k, v])
            elif isinstance(record_data, list):
                for i, item in enumerate(record_data):
                    if isinstance(item, dict):
                        for k, v in item.items():
                            writer.writerow([f"{record_name}[{i}]", k, v])
                    else:
                        writer.writerow([record_name, i, item])
            else:
                writer.writerow([record_name, "", record_data])
    print(f"✅ Saved CSV to {csv_out}")
except PermissionError:
    print(f"⚠️ Could not write {csv_out} (file may be open). Close Excel or try again.")

# --- Load parsed JSON back ---
with open(json_out, "r", encoding="utf-8") as f:
    data_dict = json.load(f)

# 🔎 Debug: show available keys and nested keys
print("\n📂 Top-level keys in parsed PDZ JSON:")
for key, val in data_dict.items():
    print("   ", key, "->", type(val))
    if isinstance(val, dict):
        print("      Nested keys:", list(val.keys())[:20])
    elif isinstance(val, list) and len(val) > 0 and isinstance(val[0], dict):
        print("      First entry keys:", list(val[0].keys())[:20])


# Peek inside Results/Elements if they exist
if "Results" in data_dict:
    print("\n📂 Sample of 'Results':")
    sample = data_dict["Results"]
    if isinstance(sample, dict):
        for k in list(sample.keys())[:10]:
            print("   ", k, ":", sample[k])
if "Elements" in data_dict:
    print("\n📂 Sample of 'Elements':")
    elements = data_dict["Elements"]
    if isinstance(elements, list):
        for e in elements[:5]:
            print("   ", e)

# --- Plot XRF Spectrum with element peak annotations ---
spectrum_info = data_dict.get("XRF Spectrum")
if spectrum_info and isinstance(spectrum_info, dict):
    spectrum = spectrum_info.get("spectrum_data", [])
    ev_per_channel = spectrum_info.get("ev_per_channel", 1)  # eV
    channel_start = spectrum_info.get("channel_start", 0)    # keV

    if spectrum:
        channels = range(len(spectrum))
        # ✅ Correct calibration: eV → keV + channel offset
        energies_keV = [(c * ev_per_channel / 1000.0) + channel_start for c in channels]

        plt.figure(figsize=(18, 6))
        plt.plot(energies_keV, spectrum, linewidth=1, label="Spectrum")

        # --- Element characteristic energies (keV) ---
        element_lines = {
            "Na": 1.04,
            "Mg": 1.25,
            "Al": 1.49,
            "S": 2.31,
            "K": 3.32,
            "Ca": 3.69,
            "Ti": 4.51,
            "V": 4.95,
            "Cr": 5.41,
            "Mn": 5.89,
            "Fe": 6.40,
            "Co": 6.92,
            "Ni": 7.47,
            "Cu": 8.04,
            "Zn": 8.63,
            "Ga": 9.24,
            "Ge": 9.87,
            "As": 10.53,
            "Se": 11.20,
            "Rb": 13.36,
            "Sr": 14.14,
            "Y": 14.92,
            "Zr": 15.74,
            "Nb": 16.57,
            "Mo": 17.43,
            "Pb": 10.55,  # Pb Lα
        }

        # --- Annotate elements only if a peak exists near the energy ---
        ymax = max(spectrum)
        threshold = 0.01 * ymax  # ignore tiny noise peaks
        tolerance_keV = 0.15     # search window around expected energy

        for elem, energy in element_lines.items():
            # Find all indices within tolerance window
            candidate_indices = [i for i, e in enumerate(energies_keV) if abs(e - energy) <= tolerance_keV]
            if not candidate_indices:
                continue

            # Find the local peak in this region
            local_vals = [spectrum[i] for i in candidate_indices]
            peak_val = max(local_vals)
            peak_idx = candidate_indices[local_vals.index(peak_val)]

            # Check threshold and annotate
            if peak_val > threshold:
                plt.text(
                    energies_keV[peak_idx],
                    spectrum[peak_idx] * 1.05,
                    elem,
                    rotation=90,
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    color="red"
                )
                print(f"🔎 {elem}: expected {energy:.2f} keV, "
                      f"found peak at {energies_keV[peak_idx]:.2f} keV "
                      f"(counts={peak_val:.0f})")

        # Keep x-axis ticks at regular intervals
        plt.xticks(range(0, 21, 2))

        plt.title("XRF Spectrum", fontsize=16)
        plt.xlim(0, 20)
        plt.xlabel("Energy (keV)", fontsize=14)
        plt.ylabel("Counts", fontsize=14)
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.tight_layout()
        plt.show()
    else:
        print("⚠️ Spectrum data not found in PDZ file.")

# --- Sherman equation placeholder ---
print("\n📐 Sherman Equation Components (from PDZ):")
if "Elements" in data_dict:
    elements = data_dict["Elements"]
    for elem in elements:
        name = elem.get("name")
        conc = elem.get("concentration")
        print(f"   {name}: C = {conc} (ppm or wt%)")
else:
    print("⚠️ No 'Elements' section found in PDZ file.")

