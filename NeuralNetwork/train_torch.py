#%% 
# import necessary libraries
import numpy as np
import torch
from torch import nn, optim
from sklearn.model_selection import train_test_split
import scipy.io
from torch.utils.data import DataLoader, TensorDataset
from stroke_model import StrokeNN

#%%
# Check if GPU is available and set device
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")
#%%
# Load and preprocess data
# Please uncomment the following to run the training on sample dataset
mat = scipy.io.loadmat('sample_training_data.mat')
X1 = mat['X1']                          # normal
X2 = mat['X2']                          # lesion
X1_clean = mat['X1_clean']              # normal
X2_clean = mat['X2_clean']              # lesion
Y1 = mat['Y1']                          # normal
Y2 = mat['Y2']                          # lesion
R11 = mat['R11']                        # normal
R12 = mat['R12']                        # lesion


# Scaling Y values to bring them in the same range
Y1[:,0] = np.log(Y1[:,0]*500)
Y1[:,1] = np.log(Y1[:,1])
Y2[:,0] = np.log(Y2[:,0]*500)
Y2[:,1] = np.log(Y2[:,1])


# Combine data
X = np.concatenate((X1, X2), axis=1)
Y = np.concatenate((Y1, Y2))
R1 = np.concatenate((R11, R12))
X = np.transpose(X[1:68, :]/R1.T)

#%%
# Train/val split
X_train, X_val, Y_train, Y_val = train_test_split(X, Y, test_size=0.3, random_state=42)

# Convert to tensors
X_train_tensor = torch.tensor(X_train, dtype=torch.float32).to(device)
Y_train_tensor = torch.tensor(Y_train, dtype=torch.float32).to(device)
X_val_tensor = torch.tensor(X_val, dtype=torch.float32).to(device)
Y_val_tensor = torch.tensor(Y_val, dtype=torch.float32).to(device)

# DataLoaders
batch_size = 64
train_loader = DataLoader(TensorDataset(X_train_tensor, Y_train_tensor), batch_size=batch_size, shuffle=True)
val_loader = DataLoader(TensorDataset(X_val_tensor, Y_val_tensor), batch_size=batch_size, shuffle=False)

input_dim = X_train.shape[1]
output_dim = Y_train.shape[1]

print(f"Training samples: {len(train_loader.dataset)}")
print(f"Validation samples: {len(val_loader.dataset)}")

#%%
# Training function
def train_model(model, train_loader, val_loader, epochs=30, lr=1e-3, patience=3):
    print(f" Starting training for {epochs} epochs...")
    model = model.to(device)
    criterion = nn.HuberLoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)
    best_val_loss = float("inf")
    patience_counter = 0

    for epoch in range(epochs):
        model.train()
        total_loss = 0

        for X_batch, Y_batch in train_loader:
            optimizer.zero_grad()
            outputs = model(X_batch)
            loss = criterion(outputs, Y_batch)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        avg_train_loss = total_loss / len(train_loader)

        model.eval()
        val_loss = 0
        with torch.no_grad():
            for X_val_batch, Y_val_batch in val_loader:
                preds = model(X_val_batch)
                val_loss += criterion(preds, Y_val_batch).item()
        avg_val_loss = val_loss / len(val_loader)

        print(f"[Epoch {epoch+1:02d}] Train Loss: {avg_train_loss:.4f} | Val Loss: {avg_val_loss:.4f}")

        # Early stopping
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            patience_counter = 0
            print(f"Validation loss improved. Saving model...")
            torch.save(model.state_dict(), 'stroke_model.pt')
        else:
            patience_counter += 1
            print(f"No improvement. Patience: {patience_counter}/{patience}")
            if patience_counter >= patience:
                print("Early stopping triggered.")
                break

    print("Training complete.")
    return model
# %%
# Initialize and train model
model = StrokeNN(input_dim=input_dim)
trained_model = train_model(model, train_loader, val_loader)

# %%
