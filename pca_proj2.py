import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
#from matplotlib.mlab import PCA as mlabPCA

path4 = "/home/prasad/Desktop/sem-3/Data Mining/Project2/iyer.txt"
path3= "/home/prasad/Desktop/sem-3/Data Mining/Project2/cho.txt"
path1 = "/home/prasad/Desktop/sem-3/Data Mining/Project2/new_dataset_1.txt"
path2 = "/home/prasad/Desktop/sem-3/Data Mining/Project2/new_dataset_2.txt"

path=path4

X = np.loadtxt(path)[:,2:]

labels1 = np.loadtxt('/home/prasad/workspace/CSE602Project2/Data/kMeans.txt')[:]
labels2 = np.loadtxt('/home/prasad/workspace/CSE602Project2/Data/Hierarchical.txt')[:]
labels3 = np.loadtxt('/home/prasad/workspace/CSE602Project2/Data/DBScan.txt')[:]
labels4 = np.loadtxt('/home/prasad/workspace/CSE602Project2/Data/mrAnalysis.txt')[:]

labels=labels3
#print labels.shape
pca = PCA(2)
new_X = pca.fit(X).transform(X)
#print new_X.shape

print('Variance Ratio: %s'
      % str(pca.explained_variance_ratio_))

plt.figure()
plt.scatter(new_X[:, 0], new_X[:, 1], c=labels, s=35.0)
plt.title('PCA')

plt.show()
