# solution_03.py

# 필요한 라이브러리 임포트
from sklearn.manifold import TSNE

# ESOL 데이터셋에서 사용할 Morgan fingerprint 배열 생성
fingerprint_array = np.array(esol['MorganFP'].tolist())  # RDKit BitVect를 numpy 배열로 변환

# t-SNE 객체 생성 (n_components=2, random_state=42 설정)
tsne = TSNE(n_components=2, random_state=42)

# t-SNE를 데이터에 적용하고 결과 저장
tsne_coordinates = tsne.fit_transform(fingerprint_array)

# 원본 데이터프레임에 tSNE 좌표 추가
esol['tSNE1'] = tsne_coordinates[:, 0]
esol['tSNE2'] = tsne_coordinates[:, 1]
