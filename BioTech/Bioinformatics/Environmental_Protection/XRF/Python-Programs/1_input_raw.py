import pandas as pd
import re

def deep_clean_line(line: str) -> str:
    """Remove nulls and collapse spaced-out characters."""
    line = line.replace("\x00", "")
    line = re.sub(r"\s+", "", line)
    return line

def normalize_token(token: str) -> str:
    """Remove whitespace inside a CSV token."""
    return re.sub(r"\s+", "", token.strip())

def extract_all_blocks(csv_file, output_file):
    blocks = []
    with open(csv_file, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    i = 0
    block_id = 0
    while i < len(lines):
        norm_line = deep_clean_line(lines[i])
        if norm_line.startswith("File#"):
            if i + 1 < len(lines):
                header = [normalize_token(tok) for tok in lines[i].replace("\x00", "").split(",")]
                row = [normalize_token(tok) for tok in lines[i+1].replace("\x00", "").split(",")]

                if len(row) < len(header):
                    row += [""] * (len(header) - len(row))

                df = pd.DataFrame([row], columns=header)
                df["__BlockID"] = block_id
                df["IncludedInAverage"] = "Yes"
                blocks.append(df)
                block_id += 1
            i += 2
        else:
            i += 1

    if not blocks:
        print("⚠️ No valid blocks found!")
        return

    final_df = pd.concat(blocks, ignore_index=True)

    # --- Reorder columns ---
    meta_end = final_df.columns.get_loc("MatchQual3") + 1 if "MatchQual3" in final_df.columns else 0
    meta_cols = final_df.columns[:meta_end].tolist()

    # Special trailing cols
    trailing_cols = [col for col in ["CalCheck", "__BlockID", "IncludedInAverage"] if col in final_df.columns]

    # Element columns = everything else
    element_cols = [c for c in final_df.columns if c not in meta_cols + trailing_cols]

    # Pair elements with their Err
    pairs = []
    used = set()
    for col in element_cols:
        if col.endswith("Err"):
            base = col[:-3]
            if base in element_cols:
                pairs.append((base, col))
                used.add(base)
                used.add(col)
    # Add unpaired (if any)
    leftovers = [c for c in element_cols if c not in used]

    # Sort alphabetically by base element name
    sorted_pairs = sorted(pairs, key=lambda x: x[0])
    element_sorted = [c for pair in sorted_pairs for c in pair] + sorted(leftovers)

    # Build final column order
    new_order = meta_cols + element_sorted + trailing_cols
    final_df = final_df[new_order]

    final_df.to_csv(output_file, index=False, encoding="utf-8")
    print(f"✅ Cleaned & reordered data saved to {output_file}")

if __name__ == "__main__":
    input_file = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Datasets\Results.csv"
    output_file = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Environmental_Protection\XRF\Datasets\Cleaned_Results.csv"

    extract_all_blocks(input_file, output_file)
