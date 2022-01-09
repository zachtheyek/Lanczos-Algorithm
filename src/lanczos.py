import numpy as np 

# Function to tri-diagonalize a matrix
def tridiag(a, b, c, k1=-1, k2=0, k3=1):
  return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

# Lanczos algorithm
def lanczos(A, v1):
  np.set_printoptions(precision=3, suppress=True)
  # Step 1
  x, y = [], []
  n = A.shape[1]
  v2, beta = 0.0, 0.0

  for i in range(n):
    # Step 2
    w_prime = np.dot(A, v1)
    conj = np.matrix.conjugate(w_prime)
    alpha = np.dot(conj, v1)
    w = w_prime - alpha * v1 - beta * v2
    # Step 3
    beta = np.linalg.norm(w)
    # Step 4
    x.append(np.linalg.norm(alpha))

    # Reset
    if i < (n-1):
        y.append(beta)
    v2 = v1
    v1 = w/beta
    
  return tridiag(y, x, y)