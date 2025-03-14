# 1. 데이터셋 읽기
data_ex = pd.read_csv('data/unknown_clusters.csv')  # 데이터셋 로드

# 2. SMILES를 Morgan fingerprints로 특징화
def featurize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)  # Morgan fingerprint 생성
    return list(fingerprint)

X_ex = data_ex['smiles'].apply(featurize_smiles)  # SMILES 컬럼을 Morgan fingerprints로 변환

# 3. 데이터를 Numpy 배열로 변환 및 표준화
X_ex = pd.DataFrame(X_ex.tolist())  # 리스트로 변환

# 4. K-means 또는 K-medoids 알고리즘으로 최적 군집 수 추정
algorithms = {
    'KMeans': KMeans(n_init=10, init='k-means++', random_state=42),
    'KMedoids': KMedoids(init='k-medoids++', metric='jaccard', random_state=42)
}

# 5. inertia와 silhouette 점수 시각화 함수 실행
plot_inertia_and_silhouette(X_ex, algorithms, min_clusters=2, max_clusters=10)

#Solution
# %load https://raw.githubusercontent.com/jfjoung/AI_For_Chemistry/main/notebooks/week5/solution_05.py