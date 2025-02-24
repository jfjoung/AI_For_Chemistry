# solution_01.py

# Importing necessary libraries
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

esol = pd.read_csv('data/esol.csv')

# Assuming `esol` DataFrame has already been loaded
# Create a 'Molecule' column containing rdkit.Mol objects from SMILES
esol['Molecule'] = esol['smiles'].apply(Chem.MolFromSmiles)

# Drop any rows where Molecule conversion failed
esol = esol.dropna(subset=['Molecule'])

# Create Morgan fingerprints (r=2, nBits=2048) from Molecule column
esol['MorganFP'] = esol['Molecule'].apply(lambda mol: AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048))
