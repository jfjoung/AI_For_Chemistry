# solution_04.py

# t-SNE 결과를 포함하여 용해도 라벨을 색상으로 구분하여 시각화
plt.figure(figsize=(8, 6))
sns.scatterplot(x=esol['tSNE1'], y=esol['tSNE2'], hue=esol['Solubility_Label'], palette='viridis', alpha=0.7)

# X, Y축 라벨 및 제목 설정
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.title('ESOL t-SNE plot')

# 범례 제목 설정
plt.legend(title='Solubility')

# 플롯 표시
plt.show()
