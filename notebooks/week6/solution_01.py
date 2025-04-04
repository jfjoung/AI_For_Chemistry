# "zinc.smi" 파일에서 학습 데이터로 사용된 SMILES 문자열을 추출합니다

zinc_smiles = []
with open('zinc.smi', 'r') as f:
    for line in f:
        smiles = line.strip()  # 줄바꿈 문자 제거
        zinc_smiles.append(smiles)

# set으로 변환해 중복 제거 + 빠른 검색 가능하도록 처리
zinc_smiles_set = set(zinc_smiles)

print(f"총 {len(zinc_smiles_set)}개의 고유한 SMILES가 학습 데이터에 포함되어 있습니다.")

#Solution
# %load https://raw.githubusercontent.com/jfjoung/AI_For_Chemistry/main/notebooks/week6/solution_01.py