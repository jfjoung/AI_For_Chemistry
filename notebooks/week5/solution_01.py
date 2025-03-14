# solution_01.py

# Importing necessary libraries
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ESOL 데이터셋 로드 
# pd.read_csv 를 사용해서 data/esol.csv에 있는 데이터를 로드해보세요.
esol = pd.read_csv('data/esol.csv')

# Chem.MolFromSmiles를 사용해 rdkit.Mol 객체를 생성하여 'Molecule' 열에 저장해보세요
# apply를 사용해서 'smiles'열에 일괄적으로 적용해보세요. 
esol['Molecule'] = esol['smiles'].apply(Chem.MolFromSmiles)

# Molecule 변환이 실패한 행 제거
# 'Molecule' 열에 있는 NaN을 제거하기 위해 dropna를 사용해보세요.
esol = esol.dropna(subset=['Molecule'])

# 'Molecule' 열을 이용하여 Morgan fingerprint (r=4, nBits=2048) 생성
# lambda mol: AllChem.GetMorganFingerprintAsBitVect(mol, radius=4, nBits=2048) 를 사용해보세요.
esol['MorganFP'] = esol['Molecule'].apply(lambda mol: AllChem.GetMorganFingerprintAsBitVect(mol, radius=4, nBits=2048))

# 결과확인
esol.head()

#Solution
# %load https://raw.githubusercontent.com/jfjoung/AI_For_Chemistry/main/notebooks/week5/solution_01.py