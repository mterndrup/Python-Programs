import cv2
import numpy as np
import os
from pathlib import Path

# === Config ===
input_dir = Path(r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\StemCells\Images")
output_dir = input_dir / "Processed"
output_dir.mkdir(exist_ok=True)

# Function to process one image
def process_image(image_path):
    img = cv2.imread(str(image_path))
    if img is None:
        print(f"⚠️ Could not read {image_path}")
        return

    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # Enhance edges using Gaussian blur + Canny
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)
    edges = cv2.Canny(blurred, 50, 150)

    # Find contours
    contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # --- First output: cells in green, edges in yellow ---
    traced = np.zeros_like(img)
    cv2.drawContours(traced, contours, -1, (0, 255, 0), 1)  # green outlines
    traced[np.where(edges > 0)] = (0, 255, 255)  # yellow edges

    # --- Second output: outlines only ---
    outlines = np.zeros_like(img)
    cv2.drawContours(outlines, contours, -1, (0, 255, 0), 1)

    # Save outputs
    base_name = image_path.stem
    cv2.imwrite(str(output_dir / f"{base_name}_traced.png"), traced)
    cv2.imwrite(str(output_dir / f"{base_name}_outlines.png"), outlines)

    print(f"✅ Processed {image_path.name}")

# === Run on all images ===
for file in input_dir.iterdir():
    if file.suffix.lower() in [".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp"]:
        process_image(file)

print(f"\nAll done! Processed images saved in: {output_dir}")
