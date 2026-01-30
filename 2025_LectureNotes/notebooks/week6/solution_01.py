from rdkit.Chem import MolToSmiles

# 생성된 Mol 객체들을 SMILES 문자열로 변환하여 리스트로 저장
generated_smiles = [MolToSmiles(mol) for mol in generated_molecules]

# 중복 제거를 위해 set으로 변환 (선택사항)
generated_smiles_set = set(generated_smiles)

print(f"총 {len(generated_smiles_set)}개의 고유한 SMILES가 생성되었습니다.")
