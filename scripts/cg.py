import numpy as np

def conjgrad(A, b, x):
    r = b - A.dot(x)
    p = r
    rsold = np.dot(r, r)
    while(True):
        Ap = A.dot(p)
        print(np.dot(p, Ap))
        alpha = rsold / np.dot(p, Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = np.dot(r, r)
        print(np.sqrt(rsnew))
        if np.sqrt(rsnew) < 1e-10:
              break
        p = r + (rsnew / rsold) * p
        rsold = rsnew
    return x

A = np.array([
    [3, 2, -1],
    [2, -2, 4],
    [-1, 4, -1],
])

b = np.array([1, -2, 0])

y = np.array([0, 0, 0])

print(conjgrad(A, b, y))
