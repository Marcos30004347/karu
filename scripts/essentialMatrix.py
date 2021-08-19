import numpy as np
from numpy.lib.twodim_base import diag

E = [
-3.410e-05,      -1.430e+03,      3.410e-05,
-1.430e+03,      1.065e-04,       1.430e+03,
-7.243e-05,      -1.430e+03,      7.243e-05,
]
E = np.reshape(E, (3,3))

U,D,V_T = np.linalg.svd(E, full_matrices=True)

e = (D[0] + D[1])/2

D[0] = e
D[1] = e
D[2] = 0

E = np.dot(np.dot(U, np.diag(D)), V_T)

U, _, V_T = np.linalg.svd(E)

W = [0,-1,0,1,0,0,0,0,1]
W = np.reshape(W, (3,3))

Z = [0,1,0,-1,0,0,0,0,0]
Z = np.reshape(Z, (3,3))

R1 = np.dot(np.dot(U, W), V_T)
R2 = np.dot(np.dot(U, W.T), V_T)

if np.linalg.det(R1) < 0:
    R1 = -R1
if np.linalg.det(R2) < 0:
    R2 = -R2

Tx = np.dot(np.dot(U, Z), U.T)

t = [Tx[2,1], Tx[0,2], Tx[1,0]]
t = np.reshape(t,(1,3))

print(t)
print(-t)

# print(R1)
# print(R2)

# print(np.dot(np.linalg.inv(R1),t.T))
# print(np.dot(np.linalg.inv(R2),t.T))
# print(np.dot(np.linalg.inv(R1),-t.T))
# print(np.dot(np.linalg.inv(R2),-t.T))

# w = [0,-1,0,1,0,0,0,0,1]

# w = np.reshape(w, (3,3))

# R1 = np.dot(np.dot(U, w), V_T)
# R2 = np.dot(np.dot(U, w.T), V_T)

# t1 = U[:, 2]
# t2 = -U[:, 2]


# pixel_length = 1
# focal_length = 1

# dx1 = 133-319.5*pixel_length
# dy1 = 75-239.5 *pixel_length
# v1  = (dx1, dy1, focal_length) - (0,0,0)

# dx2 = (124.661-319.5)*pixel_length
# dy2 = (67.6607-239.5)*pixel_length
# v2  = R1 * ( (dx2, dy2, focal_length) - (0,0,0) ) + t1

# print(np.linalg.norm(t1))
# print(np.linalg.norm(t2))
# print(t1)
# print(t2)

# C2 = np.dot(R1, t1.T)
# print(C2)
# C2 = np.dot(R1, t2.T)
# print(C2)
# C2 = np.dot(R2, t1.T)
# print(C2)
# C2 = np.dot(R2, t2.T)
# print(C2)
