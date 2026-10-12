"""Small, synthetic matrix examples. Requires NumPy; not an omics pipeline."""
import json
import numpy as np

C = np.array([[10, 20, 30], [0, 5, 10], [3, 6, 9], [8, 4, 0]])
P = np.array([[1, 0], [1, 0], [0, 1]])
group_counts = C @ P
assert group_counts.tolist() == [[30, 30], [5, 10], [9, 9], [12, 0]]

# Samples × genes. This toy log transform is not a full RNA-seq normalization.
E = np.log1p(C.T.astype(float))
E_centered = E - E.mean(axis=0, keepdims=True)
U, s, Vt = np.linalg.svd(E_centered, full_matrices=False)
scores = E_centered @ Vt[:2].T
np.testing.assert_allclose(scores, U[:, :2] * s[:2], atol=1e-12)
np.testing.assert_allclose((U * s) @ Vt, E_centered, atol=1e-12)

# One unmasked self-attention head; weights are random, not trained.
rng = np.random.default_rng(42)
X = rng.normal(size=(3, 4))
Wq = rng.normal(size=(4, 2))
Wk = rng.normal(size=(4, 2))
Wv = rng.normal(size=(4, 2))
Q, K, V = X @ Wq, X @ Wk, X @ Wv
scores_attention = Q @ K.T / np.sqrt(Q.shape[1])
z = scores_attention - scores_attention.max(axis=-1, keepdims=True)
A = np.exp(z)
A /= A.sum(axis=-1, keepdims=True)
Y = A @ V
np.testing.assert_allclose(A.sum(axis=-1), np.ones(3))
assert Y.shape == (3, 2)

if __name__ == '__main__':
    print(json.dumps({'numpy_version': np.__version__, 'group_counts':group_counts.tolist(),
        'pca_scores_shape':list(scores.shape), 'attention':A.tolist(),
        'attention_output_shape':list(Y.shape), 'status':'All example checks passed'}, indent=2))
