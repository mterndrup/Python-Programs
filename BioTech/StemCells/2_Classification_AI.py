import os
from pathlib import Path
from PIL import Image
import torch
from torch import nn, optim
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, models

# ---------------------------
# Paths
# ---------------------------
base_dir = Path(r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\StemCells")
train_dir = base_dir / "train"
inference_dir = base_dir / "Images" / "OP9"
output_file = base_dir / "classification_results.txt"
model_file = base_dir / "stem_cell_resnet18.pth"

# ---------------------------
# Create train folder if missing
# ---------------------------
train_dir.mkdir(parents=True, exist_ok=True)

# ---------------------------
# Copy all OP9 images into train folder if not already there
# ---------------------------
op9_images = [f for f in os.listdir(inference_dir) if f.lower().endswith(('.png', '.jpg', '.jpeg', '.tif', '.tiff'))]
for file in op9_images:
    src = inference_dir / file
    dst = train_dir / file
    if not dst.exists():
        dst.write_bytes(src.read_bytes())
print(f"✅ {len(op9_images)} images copied into training folder")

# ---------------------------
# Custom single-class dataset
# ---------------------------
class SingleClassDataset(Dataset):
    def __init__(self, folder, transform=None):
        self.folder = folder
        self.transform = transform
        self.img_files = [f for f in os.listdir(folder)
                          if f.lower().endswith(('.png', '.jpg', '.jpeg', '.tif', '.tiff'))]

    def __len__(self):
        return len(self.img_files)

    def __getitem__(self, idx):
        img_path = os.path.join(self.folder, self.img_files[idx])
        image = Image.open(img_path).convert("RGB")
        if self.transform:
            image = self.transform(image)
        label = 0  # all stem_cell
        return image, label

# ---------------------------
# Device
# ---------------------------
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# ---------------------------
# Data transformations with augmentation
# ---------------------------
train_transform = transforms.Compose([
    transforms.Resize((224,224)),
    transforms.RandomHorizontalFlip(),
    transforms.RandomRotation(20),
    transforms.ColorJitter(brightness=0.2, contrast=0.2),
    transforms.ToTensor(),
    transforms.Normalize([0.485,0.456,0.406],[0.229,0.224,0.225])
])

# ---------------------------
# Dataset and loader
# ---------------------------
train_dataset = SingleClassDataset(train_dir, transform=train_transform)
train_loader = DataLoader(train_dataset, batch_size=8, shuffle=True)

# ---------------------------
# Model setup
# ---------------------------
model = models.resnet18(weights=models.ResNet18_Weights.DEFAULT)
num_features = model.fc.in_features
model.fc = nn.Linear(num_features, 1)  # single output
model = model.to(device)

criterion = nn.BCEWithLogitsLoss()
optimizer = optim.Adam(model.parameters(), lr=1e-4)

# ---------------------------
# Training loop
# ---------------------------
num_epochs = 10  # train longer for better results

for epoch in range(num_epochs):
    model.train()
    running_loss = 0.0
    for inputs, labels in train_loader:
        inputs = inputs.to(device)
        labels = labels.float().unsqueeze(1).to(device)

        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        running_loss += loss.item() * inputs.size(0)

    epoch_loss = running_loss / len(train_dataset)
    print(f"Epoch {epoch+1}/{num_epochs} - Loss: {epoch_loss:.4f}")

# ---------------------------
# Save model
# ---------------------------
torch.save(model.state_dict(), model_file)
print(f"✅ Model saved to {model_file}")

# ---------------------------
# Inference on OP9 images
# ---------------------------
model.eval()
results = []

val_transform = transforms.Compose([
    transforms.Resize((224,224)),
    transforms.ToTensor(),
    transforms.Normalize([0.485,0.456,0.406],[0.229,0.224,0.225])
])

for file in op9_images:
    img_path = inference_dir / file
    image = Image.open(img_path).convert("RGB")
    input_tensor = val_transform(image).unsqueeze(0).to(device)
    with torch.no_grad():
        output = model(input_tensor)
        pred = torch.sigmoid(output)
        label = "stem_cell" if pred.item() > 0.5 else "unknown"
    line = f"{file}: {label}"
    print(line)
    results.append(line)

# Save inference results
with open(output_file, "w") as f:
    f.write("\n".join(results))
print(f"\n✅ Classification complete. Results saved to {output_file}")
