__author__ = 'Michael'
import numpy as np
from sklearn.decomposition import PCA



def main():
    print("Hello World!")
    X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    pca = PCA(n_components=1)
    pca.fit(X)
    print(pca.explained_variance_ratio_)
    print(pca.transform(X))


main()