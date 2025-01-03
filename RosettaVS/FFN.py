import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, log_loss
from sklearn.preprocessing import StandardScaler
import torch.nn.functional as F

def replace_top_percent(arr, threshold):
    sorted_arr = np.sort(arr)
    threshold_value = sorted_arr[int(len(arr) * (1 - threshold))]
    result = np.where(arr >= threshold_value, 1, 0)
    return result

def generate_morgan_fingerprint(smiles, radius=2, n_bits=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits)
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)

class MolecularDataset(Dataset):
    def __init__(self, smiles_list, labels, radius=2, n_bits=1024):
        self.smiles_list = smiles_list
        self.labels = labels
        self.radius = radius
        self.n_bits = n_bits

    def __len__(self):
        return len(self.smiles_list)

    def __getitem__(self, idx):
        smiles = self.smiles_list[idx]
        label = self.labels[idx]
        features = generate_morgan_fingerprint(smiles, self.radius, self.n_bits)
        return torch.tensor(features, dtype=torch.float32), torch.tensor(label, dtype=torch.long)

class MolecularModel(nn.Module):
    def __init__(self, input_size=1024):
        super(MolecularModel, self).__init__()
        self.fc1 = nn.Linear(input_size, 512)
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, 2) 

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)
        return x


def train_model(model, train_loader, val_loader, epochs=50, lr=0.001, patience=10):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    
    criterion = nn.CrossEntropyLoss() 
    optimizer = optim.Adam(model.parameters(), lr=lr)
    
    best_loss = np.inf
    patience_counter = 0
    best_model_state = None

    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        correct_preds = 0
        total_preds = 0

        for inputs, labels in train_loader:
            inputs, labels = inputs.to(device), labels.to(device)

            optimizer.zero_grad()

            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()

            running_loss += loss.item()

            _, predicted = torch.max(outputs.data, 1)
            total_preds += labels.size(0)
            correct_preds += (predicted == labels).sum().item()

        avg_train_loss = running_loss / len(train_loader)
        train_accuracy = correct_preds / total_preds

        val_loss, val_accuracy = evaluate_model(model, val_loader, device)

        print(f"Epoch [{epoch+1}/{epochs}] | Train Loss: {avg_train_loss:.4f} | Train Accuracy: {train_accuracy:.4f} | Validation Loss: {val_loss:.4f} | Validation Accuracy: {val_accuracy:.4f}")

        # Early stopping
        if val_loss < best_loss:
            best_loss = val_loss
            best_model_state = model.state_dict()
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"Early stopping triggered. Stopping at epoch {epoch+1}.")
                break

    model.load_state_dict(best_model_state)
    return model

def evaluate_model(model, val_loader, device):
    model.eval()
    val_loss = 0.0
    correct_preds = 0
    total_preds = 0
    all_labels = []
    all_preds = []

    with torch.no_grad():
        for inputs, labels in val_loader:
            inputs, labels = inputs.to(device), labels.to(device)

            outputs = model(inputs)
            criterion = nn.CrossEntropyLoss()
            loss = criterion(outputs, labels)
            val_loss += loss.item()

            _, predicted = torch.max(outputs.data, 1)
            total_preds += labels.size(0)
            correct_preds += (predicted == labels).sum().item()

            all_labels.extend(labels.cpu().numpy())
            all_preds.extend(predicted.cpu().numpy())

    avg_val_loss = val_loss / len(val_loader)
    val_accuracy = correct_preds / total_preds

    return avg_val_loss, val_accuracy



def get_model(df_train, df_test, dG_cutoff):
    smileses_train = df_train['Smiles'].values
    dGs_train = df_train['dG'].values
    smileses_test = df_test['Smiles'].values
    dGs_test = df_test['dG'].values
    labels_train = replace_top_percent(dGs_train, dG_cutoff)
    labels_test = replace_top_percent(dGs_test, dG_cutoff)

    train_dataset = MolecularDataset(smileses_train, labels_train)
    val_dataset = MolecularDataset(smileses_test, labels_test)

    train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=64, shuffle=False)

    model = MolecularModel(input_size=1024)
    best_model = train_model(model, train_loader, val_loader, epochs=100, lr=0.001, patience=10)
    return best_model



def predict(model, df_screen, device):
    model.eval()
    all_probs = [] 
    smileses_train = df_screen['Smiles'].values
    pred_dataset = MolecularDataset(smileses_train, [0.5 for i in smileses_train])
    pred_loader = DataLoader(pred_dataset, batch_size=64, shuffle=False)  # shuffle=False for evaluation
    
    with torch.no_grad():
        for inputs, labels in pred_loader:
            inputs = inputs.to(device)
            outputs = model(inputs)
            probs = F.softmax(outputs, dim=1)
            prob_class_1 = probs[:, 1].cpu().numpy()  
            all_probs.extend(prob_class_1)
    df_screen['Probability'] = all_probs
    return df_screen



