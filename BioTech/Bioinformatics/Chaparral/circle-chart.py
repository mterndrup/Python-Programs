import matplotlib.pyplot as plt
import numpy as np

# --- Inner Circle: Wildfire Succession Phases ---
inner_labels = [
    "Acute Burn\nPhase",
    "Soil\nStabilization",
    "Transitional\nBuild-up",
    "Pioneer\nColonization"
]
inner_sizes = [10, 35, 35, 20]
inner_colors = ['#EDAEAE', '#EDAEAE', '#EDAEAE', '#EDAEAE']

# --- Outer Circle: Functional Groups ---
outer_labels = [
    "Heat-Activated Germinators",
    "Mycorrhizal\nRe-Establishers",
    "Plant Debris Composers",
    "PyOM Degraders",
    "Rapid\nColonizers"
]
outer_sizes = [10, 20, 25, 25, 20]
outer_colors = ['#2D8F65', '#2D8F65', '#2D8F65', '#2D8F65', '#2D8F65']

fig, ax = plt.subplots(figsize=(8, 8))

# Draw outer pie WITHOUT labels (we'll add them manually)
wedges_outer, _ = ax.pie(
    outer_sizes,
    radius=1.0,
    labels=[None] * len(outer_labels),
    colors=outer_colors,
    startangle=90,
    wedgeprops=dict(width=0.5, edgecolor='white')
)

# Draw inner pie WITHOUT labels
wedges_inner, _ = ax.pie(
    inner_sizes,
    radius=0.5,
    labels=[None] * len(inner_labels),
    colors=inner_colors,
    startangle=90,
    wedgeprops=dict(width=0.5, edgecolor='white')
)

# --- Calculate angles for outer wedges to position labels ---
total_outer = sum(outer_sizes)
angles_outer = np.cumsum([0] + outer_sizes) / total_outer * 360
start_angle = 90

for i, label in enumerate(outer_labels):
    angle = (angles_outer[i] + angles_outer[i + 1]) / 2
    angle = (start_angle - angle) % 360
    theta = np.deg2rad(angle)

    # Base radius for outer labels (slightly outside the outer circle)
    r = 1.0 + 0.15

    # Default position
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # Custom adjustments for each outer label
    if label == "Heat-Activated\nGerminators":
        x -= 0.1
        y += 0.05
    elif label == "Mycorrhizal\nRe-Establishers":
        x += -1.8
        y -=  0.15
    elif label == "Plant Debris Composers":
        y -= 0.15
    elif label == "PyOM Degraders":
        x -= 0.05
        y -= 0.2
    elif label == "Rapid\nColonizers":
        x += 0.2

    ax.text(x, y, label, ha='center', va='center', fontsize=9, fontweight='bold')

# --- Calculate angles for inner wedges to position labels ---
total_inner = sum(inner_sizes)
angles_inner = np.cumsum([0] + inner_sizes) / total_inner * 360

for i, label in enumerate(inner_labels):
    angle = (angles_inner[i] + angles_inner[i + 1]) / 2
    angle = (start_angle - angle) % 360
    theta = np.deg2rad(angle)

    # Base radius for inner labels
    r = 0.5 - 0.2  # radius 0.5 minus width/2 of inner ring (0.5 width)

    # Default position
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # Custom adjustments for each inner label
    if label == "Acute Burn\nPhase":
        x -= 0.15
    elif label == "Pioneer\nColonization":
        x += 0.4
    elif label == "Soil\nStabilization":
        y -= 0.15
        x -= 0.1
    elif label == "Transitional\nBuild-up":
        y = 0.025
        x -= 0.05

    ax.text(x, y, label, ha='center', va='center', fontsize=9, fontweight='bold')

ax.set(aspect='equal')
plt.title("Wildfire Succession Phases and Functional Groups")
plt.show()
