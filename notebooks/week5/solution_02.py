# solution_02.py

# 용해도 값을 기준으로 카테고리 라벨을 생성하는 함수
def solubility_class(log_sol):
    '''용해도 값에 따라 해당하는 라벨을 반환하는 함수'''
    if log_sol < -5:
        return 'Low'
    elif log_sol > -1:
        return 'High'
    else:
        return 'Medium'

### YOUR CODE ####
# 'log solubility (mol/L)' 값을 이용해 'Solubility_Label' 열에 라벨 추가
esol['Solubility_Label'] = esol['log solubility (mol/L)'].apply(solubility_class)

# 새로운 라벨을 포함한 PCA 플롯 생성
plt.figure(figsize=(8, 6))
sns.scatterplot(x=esol['PC1'], y=esol['PC2'], hue=esol['Solubility_Label'], palette='viridis', alpha=0.7)

# X, Y축 라벨 및 제목 설정
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('ESOL PCA plot with solubility label')

# 범례 제목 설정
plt.legend(title='Solubility')

# 플롯 표시
plt.show()


#Solution
# %load https://raw.githubusercontent.com/jfjoung/AI_For_Chemistry/main/notebooks/week5/solution_02.py


