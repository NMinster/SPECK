import numpy as np
from scipy.linalg import svd
from sklearn.cluster import KMeans

#There is no rle and Ckmeans libraries in python so we need to define a function
#Keep in mind this can produce inconsistent results if there is no random seed set
def rle(inarray):
    ia = np.asarray(inarray)
    n = len(ia)
    if n == 0: 
        return None, None, None, None

    y = np.array(ia[1:] != ia[:-1])           # pairwise unequal 
    i = np.append(np.where(y), n - 1)         # must include last element
    lengths = np.diff(np.append(-1, i))       # run lengths
    values = ia[i]                            # run values
    
    return values, lengths

# input matrix
Xm_n = np.random.rand(100, 100)  # just an example

# Step 1: Normalize the input
Xm_n_normalized = Xm_n / np.linalg.norm(Xm_n)

# Step 2: SVD
U, s, Vh = svd(Xm_n_normalized)

# Select rank-100 SVD
mu = U[:, :100]
sigma = np.diag(s[:100])
v = Vh[:100, :]

# Step 3: Compute PCs stdev
rsdev = np.sqrt(sigma / (len(mu) - 1))

# Step 4: Calculate rate of change
rdiff = np.abs(np.diff(rsdev, axis=0))

# Step 5 and 6: Run length encoding
rval, rlen = rle(rdiff)

# Step 7: Rank-k
k = np.min(np.where((rval >= 0.01) & (rlen >= 2)))

# Step 8: Matrix reconstruction
Xm_n = mu[:, :k] @ sigma[:k, :k] @ v[:k, :]

# Step 9 to 16: Clustered thresholding
n = Xm_n.shape[1]
for v in range(n):
    X = Xm_n[v, :].reshape(-1, 1)
    kmeans = KMeans(n_clusters=4, random_state=0).fit(X)
    num = kmeans.labels_
    val = kmeans.cluster_centers_
    
    if len(set(num)) > 1:
        ind = np.argmin(val)
        Xm_n[v, num == ind] = 0
