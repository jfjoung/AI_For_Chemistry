# 세미콜론 구분자로 된 DFT 특성 CSV 파일 불러오기
df_dft = pd.read_csv('data/Dimer_LKB_P.csv', sep=';')

# 데이터프레임의 기본 정보 확인
df_dft.head()