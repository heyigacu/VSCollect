import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import pairwise_distances
from sklearn.cluster import AgglomerativeClustering

# 1. Load SMILES data from a CSV file
def load_smiles_from_csv(file_path):
    df = pd.read_csv(file_path)
    smiles_list = df['Smiles'].dropna().tolist()  # Ensure no NaN values
    return smiles_list

# 2. Convert SMILES to RDKit molecule objects
def smiles_to_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol

# 3. Calculate ECFP (Morgan Fingerprints) for molecules
def calculate_fingerprints(molecules):
    fingerprints = []
    for mol in molecules:
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fingerprints.append(fp)
    return fingerprints

# 4. Calculate Tanimoto similarity matrix
def calculate_similarity(fingerprints):
    bit_vectors = [np.array(fp) for fp in fingerprints]
    similarity_matrix = pairwise_distances(bit_vectors, metric="jaccard")
    similarity_matrix = 1 - similarity_matrix  # Convert to similarity (1 - distance)
    return similarity_matrix

# 5. Perform clustering using hierarchical clustering
def cluster_molecules(similarity_matrix, threshold=0.6):
    clustering = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='complete', distance_threshold=1-threshold)
    clustering.fit(1 - similarity_matrix)  # 1 - similarity gives distance matrix
    return clustering.labels_

# 6. Select centroids (smallest molecule in each cluster)
def select_centroids(smiles_list, labels):
    centroids = []
    for label in np.unique(labels):
        cluster_indices = np.where(labels == label)[0]
        cluster_smiles = [smiles_list[i] for i in cluster_indices]
        # Select the smallest molecule by SMILES string length (or other criteria)
        centroid = min(cluster_smiles, key=len)
        centroids.append(centroid)
    return centroids

# 7. Main process to create the druglike-centroid library
def create_druglike_centroid_library(file_path, threshold=0.6):
    smiles_list = load_smiles_from_csv(file_path)
    molecules = [smiles_to_molecule(smiles) for smiles in smiles_list]
    fingerprints = calculate_fingerprints(molecules)
    similarity_matrix = calculate_similarity(fingerprints)
    labels = cluster_molecules(similarity_matrix, threshold)
    centroids = select_centroids(smiles_list, labels)
    return centroids

# 8. Save centroids to a new CSV file
def save_centroids_to_csv(centroids, output_file):
    df_centroids = pd.DataFrame(centroids, columns=['Smiles'])
    df_centroids.to_csv(output_file, index=False)

# Example usage
file_path = 'zinc_library.csv'
output_file = 'druglike_centroid_library.csv'
centroid_library = create_druglike_centroid_library(file_path, threshold=0.6)
save_centroids_to_csv(centroid_library, output_file)
print(f"Centroids saved to {output_file}")
