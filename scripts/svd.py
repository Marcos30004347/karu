import numpy as np
from sympy import *
import math

A = [
2.207e+00, 7.848e-01, 1.529e+00, 1.044e+00, 3.713e-01, 7.235e-01, 1.444e+00, 5.133e-01, 1.000e+00,
2.781e+00, -1.289e+00,-1.737e+00,-8.225e-01,3.812e-01, 5.138e-01, -1.601e+00,7.419e-01, 1.000e+00,
5.524e-01, -4.081e-01,-6.306e-01,-4.821e-01,3.562e-01, 5.504e-01, -8.759e-01,6.472e-01, 1.000e+00,
8.665e-01, 4.610e-01, 8.388e-01, 7.075e-01, 3.764e-01, 6.849e-01, 1.033e+00, 5.496e-01, 1.000e+00,
2.207e+00, -9.372e-01,1.529e+00, -8.923e-01,3.789e-01, -6.181e-01,1.444e+00, -6.130e-01,1.000e+00,
2.781e+00, 1.065e+00, -1.737e+00,9.895e-01, 3.789e-01, -6.181e-01,-1.601e+00,-6.130e-01,1.000e+00,
5.524e-01, 3.866e-01, -6.306e-01,5.414e-01, 3.789e-01, -6.181e-01,-8.759e-01,-6.130e-01,1.000e+00,
8.665e-01, -5.142e-01,8.388e-01, -6.386e-01,3.789e-01, -6.181e-01,1.033e+00, -6.130e-01,1.000e+00,
]

A = np.reshape(A, (8,9))

print(np.linalg.svd(A)[0])
print()
print(np.linalg.svd(A)[1])
print()
print(np.linalg.svd(A)[2])
exit()
# A = [
#  [-0.00310695, -0.0025646,    2.96584],
#  [-0.028094,   -0.00771621,   56.3813],
#  [13.1905,     -29.2007,     -9999.79],
# ]

# A = Matrix([
# [9.091e-02,       -9.091e-02,      -2.727e-01,      -9.091e-02,      9.091e-02,       2.727e-01,       -3.333e-01,      3.333e-01,       1.000e+00],
# [-7.438e-02,      7.438e-02,       2.727e-01,       -7.438e-02,      7.438e-02,       2.727e-01,       -2.727e-01,      2.727e-01,       1.000e+00],
# [9.091e-02,       9.091e-02,       3.333e-01,       9.091e-02,       9.091e-02,       3.333e-01,       2.727e-01,       2.727e-01,       1.000e+00],
# [-1.111e-01,      -1.111e-01,      -3.333e-01,      1.111e-01,       1.111e-01,       3.333e-01,       3.333e-01,       3.333e-01,       1.000e+00],
# [9.091e-02,       -2.168e-08,      -2.727e-01,      2.679e-08,       -6.389e-15,      -8.038e-08,      -3.333e-01,      7.948e-08,       1.000e+00],
# [-7.438e-02,      2.128e-08,       2.727e-01,       2.192e-08,       -6.273e-15,      -8.038e-08,      -2.727e-01,      7.803e-08,       1.000e+00],
# [9.091e-02,       2.168e-08,       3.333e-01,       -2.192e-08,      -5.227e-15,      -8.038e-08,      2.727e-01,       6.503e-08,       1.000e+00],
# [-1.111e-01,      -2.119e-08,      -3.333e-01,      -2.679e-08,      -5.111e-15,      -8.038e-08,      3.333e-01,       6.358e-08,       1.000e+00],
# ])

# # print(A[0])
# # print(A[1])
# # print(A[2])
# # print(A[3])
# # print(A[4])
# # print(A[5])
# # print(A[6])
# # print(A[7])
# M, v = A.rref()
# M = np.reshape(M, (8,9))

# print(M)





# def rref(B, tol=1e-8, debug=False):
#   A = B.copy()

#   rows, cols = A.shape

#   r = 0

#   pivots_pos = []

#   row_exchanges = np.arange(rows)

#   for c in range(cols):

#     ## Find the pivot row:
#     pivot = np.argmax (np.abs (A[r:rows,c])) + r
#     m = np.abs(A[pivot, c])

#     if debug: 
#         print("Found pivot", m, "in row", pivot)

#     if m <= tol:
#       ## Skip column c, making sure the approximately zero terms are
#       ## actually zero.
#       A[r:rows, c] = np.zeros(rows-r)
      
#       if debug: 
#         print("All elements at and below (", r, ",", c, ") are zero.. moving on..")
    
#     else:
#       ## keep track of bound variables
#       pivots_pos.append((r,c))

#       if pivot != r:
#         ## Swap current row and pivot row
#         A[[pivot, r], c:cols] = A[[r, pivot], c:cols]
#         row_exchanges[[pivot,r]] = row_exchanges[[r,pivot]]
        
#         if debug: 
#             print("Swap row", r, "with row", pivot, "Now:")
#             print(A)

#       ## Normalize pivot row
#       A[r, c:cols] = A[r, c:cols] / A[r, c];

#       ## Eliminate the current column
#       v = A[r, c:cols]
#       ## Above (before row r):
#       if r > 0:
#         ridx_above = np.arange(r)
#         A[ridx_above, c:cols] = A[ridx_above, c:cols] - np.outer(v, A[ridx_above, c]).T
#         if debug:
#             print("Elimination above performed:");
#             print(A)
#       ## Below (after row r):
#       if r < rows-1:
#         ridx_below = np.arange(r+1,rows)
#         A[ridx_below, c:cols] = A[ridx_below, c:cols] - np.outer(v, A[ridx_below, c]).T
#         if debug:
#             print("Elimination below performed:")
#             print(A)
#       r += 1
#     ## Check if done
#     if r == rows:
#       break;
#   return (A, pivots_pos, row_exchanges)




# A = np.reshape(A, (8, 9))

# print(rref(A, debug=True))


# print(rref(A, debug=True))
# v = [-M[0][8], -M[1][8], -M[2][8], -M[3][8], -M[4][8], -M[5][8], -M[6][8], -M[7][8], 1]
# print(v)

# v = np.reshape(v, (1,9))

# print(np.dot(A, v.T))


Q = [
	-6.201e-14,      -2.601e-06,      1.860e-10,
	3.910e-06,       -9.765e-05,      7.802e-03,
	5.939e-10,       1.173e-02,       -2.118e-07,
]

Q = np.reshape(Q, (3,3))

# print(np.linalg.det(A))
U,D,W = np.linalg.svd(Q, full_matrices=True)

print(D)
D[2] = 0

Q = np.matmul(np.matmul(U, np.diag(D)), W)

U,D,W_T = np.linalg.svd(Q, full_matrices=True)

W = W_T.T

if np.linalg.det(W) < 0:
	W = -W

if np.linalg.det(U) < 0:
	U = -U

print(D)

r = D[0]
s = D[1]

E = [0, 1, 0,-1, 0, 0,0, 0, 1]

E = np.reshape(E, (3,3))

V = np.matmul(W,E)

print(np.linalg.det(U))
print(np.linalg.det(W))
print(np.linalg.det(V))

M1 = Matrix([
    [ Symbol("U[0][0]")*Symbol("V[0][2]"), Symbol("U[0][1]")*Symbol("V[0][2]"), Symbol("U[0][2]")*Symbol("V[0][2]"), Symbol("r")*Symbol("U[0][0]")*Symbol("V[0][0]") + Symbol("s")*Symbol("U[0][1]")*Symbol("V[0][1]") ],
    [ Symbol("U[0][0]")*Symbol("V[1][2]"), Symbol("U[0][1]")*Symbol("V[1][2]"), Symbol("U[0][2]")*Symbol("V[1][2]"), Symbol("r")*Symbol("U[0][0]")*Symbol("V[1][0]") + Symbol("s")*Symbol("U[0][1]")*Symbol("V[1][1]") ],
    [ Symbol("U[1][0]")*Symbol("V[0][2]"), Symbol("U[1][1]")*Symbol("V[0][2]"), Symbol("U[1][2]")*Symbol("V[0][2]"), Symbol("r")*Symbol("U[1][0]")*Symbol("V[0][0]") + Symbol("s")*Symbol("U[1][1]")*Symbol("V[0][1]") ],
    [ Symbol("U[1][0]")*Symbol("V[1][2]"), Symbol("U[1][1]")*Symbol("V[1][2]"), Symbol("U[1][2]")*Symbol("V[1][2]"), Symbol("r")*Symbol("U[1][0]")*Symbol("V[1][0]") + Symbol("s")*Symbol("U[1][1]")*Symbol("V[1][1]") ],
])

Mx = Matrix([
    [-Symbol("s")*Symbol("U[0][2]")*Symbol("V[0][0]"), -Symbol("r")*Symbol("U[0][2]")*Symbol("V[0][1]"), Symbol("r")*Symbol("U[0][1]")*Symbol("V[0][1]") + Symbol("s")*Symbol("U[0][0]")*Symbol("V[0][0]"), Symbol("r")*Symbol("s")*Symbol("U[0][2]")*Symbol("V[0][2]")],
    [-Symbol("s")*Symbol("U[0][2]")*Symbol("V[1][0]"), -Symbol("r")*Symbol("U[0][2]")*Symbol("V[1][1]"), Symbol("r")*Symbol("U[0][1]")*Symbol("V[1][1]") + Symbol("s")*Symbol("U[0][0]")*Symbol("V[1][0]"), Symbol("r")*Symbol("s")*Symbol("U[0][2]")*Symbol("V[1][2]")],
    [-Symbol("s")*Symbol("U[1][2]")*Symbol("V[0][0]"), -Symbol("r")*Symbol("U[1][2]")*Symbol("V[0][1]"), Symbol("r")*Symbol("U[1][1]")*Symbol("V[0][1]") + Symbol("s")*Symbol("U[1][0]")*Symbol("V[0][0]"), Symbol("r")*Symbol("s")*Symbol("U[1][2]")*Symbol("V[0][2]")],
    [-Symbol("s")*Symbol("U[1][2]")*Symbol("V[1][0]"), -Symbol("r")*Symbol("U[1][2]")*Symbol("V[1][1]"), Symbol("r")*Symbol("U[1][1]")*Symbol("V[1][1]") + Symbol("s")*Symbol("U[1][0]")*Symbol("V[1][0]"), Symbol("r")*Symbol("s")*Symbol("U[1][2]")*Symbol("V[1][2]")],
])

M_ = Matrix(M1 - Symbol("x")*Mx)

deter = det(M_)

px = Poly(deter, Symbol("x"))

coeffs = px.all_coeffs()

def getPolyRoot(U, V, r, s):
	px_ = px.subs(Symbol("U[0][0]"), U[0][0])
	px_ = px_.subs(Symbol("U[0][1]"), U[0][1])
	px_ = px_.subs(Symbol("U[0][2]"), U[0][2])
	px_ = px_.subs(Symbol("U[1][0]"), U[1][0])
	px_ = px_.subs(Symbol("U[1][1]"), U[1][1])
	px_ = px_.subs(Symbol("U[1][2]"), U[1][2])
	px_ = px_.subs(Symbol("U[2][0]"), U[2][0])
	px_ = px_.subs(Symbol("U[2][1]"), U[2][1])
	px_ = px_.subs(Symbol("U[2][2]"), U[2][2])
	px_ = px_.subs(Symbol("V[0][0]"), V[0][0])
	px_ = px_.subs(Symbol("V[0][1]"), V[0][1])
	px_ = px_.subs(Symbol("V[0][2]"), V[0][2])
	px_ = px_.subs(Symbol("V[1][0]"), V[1][0])
	px_ = px_.subs(Symbol("V[1][1]"), V[1][1])
	px_ = px_.subs(Symbol("V[1][2]"), V[1][2])
	px_ = px_.subs(Symbol("V[2][0]"), V[2][0])
	px_ = px_.subs(Symbol("V[2][1]"), V[2][1])
	px_ = px_.subs(Symbol("V[2][2]"), V[2][2])
	px_ = px_.subs(Symbol("r"), r)
	px_ = px_.subs(Symbol("s"), s)
	px_ = Poly(px_, Symbol("x"))
	print(roots(px_))
	c = px_.all_coeffs()

	# _M1 = M1.subs(Symbol("U[0][0]"), U[0][0])
	# _M1 = _M1.subs(Symbol("U[0][1]"), U[0][1])
	# _M1 = _M1.subs(Symbol("U[0][2]"), U[0][2])
	# _M1 = _M1.subs(Symbol("U[1][0]"), U[1][0])
	# _M1 = _M1.subs(Symbol("U[1][1]"), U[1][1])
	# _M1 = _M1.subs(Symbol("U[1][2]"), U[1][2])
	# _M1 = _M1.subs(Symbol("U[2][0]"), U[2][0])
	# _M1 = _M1.subs(Symbol("U[2][1]"), U[2][1])
	# _M1 = _M1.subs(Symbol("U[2][2]"), U[2][2])
	# _M1 = _M1.subs(Symbol("V[0][0]"), V[0][0])
	# _M1 = _M1.subs(Symbol("V[0][1]"), V[0][1])
	# _M1 = _M1.subs(Symbol("V[0][2]"), V[0][2])
	# _M1 = _M1.subs(Symbol("V[1][0]"), V[1][0])
	# _M1 = _M1.subs(Symbol("V[1][1]"), V[1][1])
	# _M1 = _M1.subs(Symbol("V[1][2]"), V[1][2])
	# _M1 = _M1.subs(Symbol("V[2][0]"), V[2][0])
	# _M1 = _M1.subs(Symbol("V[2][1]"), V[2][1])
	# _M1 = _M1.subs(Symbol("V[2][2]"), V[2][2])
	# _M1 = _M1.subs(Symbol("r"), r)
	# _M1 = _M1.subs(Symbol("s"), s)

	# _Mx = Mx.subs(Symbol("U[0][0]"), U[0][0])
	# _Mx = _Mx.subs(Symbol("U[0][1]"), U[0][1])
	# _Mx = _Mx.subs(Symbol("U[0][2]"), U[0][2])
	# _Mx = _Mx.subs(Symbol("U[1][0]"), U[1][0])
	# _Mx = _Mx.subs(Symbol("U[1][1]"), U[1][1])
	# _Mx = _Mx.subs(Symbol("U[1][2]"), U[1][2])
	# _Mx = _Mx.subs(Symbol("U[2][0]"), U[2][0])
	# _Mx = _Mx.subs(Symbol("U[2][1]"), U[2][1])
	# _Mx = _Mx.subs(Symbol("U[2][2]"), U[2][2])
	# _Mx = _Mx.subs(Symbol("V[0][0]"), V[0][0])
	# _Mx = _Mx.subs(Symbol("V[0][1]"), V[0][1])
	# _Mx = _Mx.subs(Symbol("V[0][2]"), V[0][2])
	# _Mx = _Mx.subs(Symbol("V[1][0]"), V[1][0])
	# _Mx = _Mx.subs(Symbol("V[1][1]"), V[1][1])
	# _Mx = _Mx.subs(Symbol("V[1][2]"), V[1][2])
	# _Mx = _Mx.subs(Symbol("V[2][0]"), V[2][0])
	# _Mx = _Mx.subs(Symbol("V[2][1]"), V[2][1])
	# _Mx = _Mx.subs(Symbol("V[2][2]"), V[2][2])
	# _Mx = _Mx.subs(Symbol("r"), r)
	# _Mx = _Mx.subs(Symbol("s"), s)

	# M = Matrix(_M1 - Symbol("x")*_Mx)
	# # print(M)
	# dt = det(M)
	# n,d = fraction(dt)
	# print(n)
	# print(d)
	# n = Poly(n, Symbol("x"))
	# d = Poly(d, Symbol("x"))
	# q, r = pdiv(n,d)
	# print(q)
	# print(r)
	# c = q.all_coeffs()

	# print(c)
	# print(roots(px_))

	return math.sqrt(-c[2]/c[0])


def solveSystem(U, V, r, s, x):
	M = Matrix(M1 - x*Mx)
	M = M.subs(Symbol("U[0][0]"), U[0][0])
	M = M.subs(Symbol("U[0][1]"), U[0][1])
	M = M.subs(Symbol("U[0][2]"), U[0][2])
	M = M.subs(Symbol("U[1][0]"), U[1][0])
	M = M.subs(Symbol("U[1][1]"), U[1][1])
	M = M.subs(Symbol("U[1][2]"), U[1][2])
	M = M.subs(Symbol("U[2][0]"), U[2][0])
	M = M.subs(Symbol("U[2][1]"), U[2][1])
	M = M.subs(Symbol("U[2][2]"), U[2][2])
	M = M.subs(Symbol("V[0][0]"), V[0][0])
	M = M.subs(Symbol("V[0][1]"), V[0][1])
	M = M.subs(Symbol("V[0][2]"), V[0][2])
	M = M.subs(Symbol("V[1][0]"), V[1][0])
	M = M.subs(Symbol("V[1][1]"), V[1][1])
	M = M.subs(Symbol("V[1][2]"), V[1][2])
	M = M.subs(Symbol("V[2][0]"), V[2][0])
	M = M.subs(Symbol("V[2][1]"), V[2][1])
	M = M.subs(Symbol("V[2][2]"), V[2][2])
	M = M.subs(Symbol("r"), r)
	M = M.subs(Symbol("s"), s)
	M = M.subs(Symbol("x"), x)

	M = Matrix(M).row_join(Matrix([[0],[0],[0],[0]]))

	# M = np.reshape(M, (4,4))
	# M = M.astype(float)
	# M = np.reshape(M, (4,4))
	print(M.rref())
	ns = M.nullspace()
	print(ns)
	M = np.array(M).astype(np.float64)
	_,Z,C = np.linalg.svd(M)
	print(C)
	return (C[3][0], C[3][1], C[3][2], C[3][3])

x = getPolyRoot(U,V,r,s)
print(x)
a, b, y, _ = solveSystem(U, V, r, s, x)

X = [r, 0, a, 0, s, b, 0, 0, y]
X = np.reshape(X, (3,3))

Xstar = [s*y, 0, 0, 0, r*y, 0, -s*a, -r*b, r*s]
Xstar = np.reshape(Xstar, (3,3))

f = np.dot(U, X)
f = np.dot(f, V.T)

g = np.dot(U, Xstar)
g = np.dot(g, V.T)

print(f)
print(g)

k2 = sqrt(x*g[2][0]/f[2][0])
k1 = sqrt(f[0][2]/(x*g[0][2]))

print(k1)
print(k2)

K1 = [
    1, 0, 0,
    0, 1, 0,
    0, 0, k1
]
K2 = [
    1, 0, 0,
    0, 1, 0,
    0, 0, k2
]

K1 = np.reshape(K1,(3,3))
K2 = np.reshape(K2,(3,3))

A = np.dot(np.dot(K2,Q), K1)
A = A.astype(float)
A = np.reshape(A, (3,3))

U_,D_,V_ = np.linalg.svd(A, full_matrices=True)

t1 = np.zeros((1,3))
P1 = np.concatenate((K1, t1.T), axis=1)


R21 = np.dot(np.dot(np.dot(K2, U_), E), V_)
R22 = np.dot(np.dot(np.dot(K2, U_), E.T),V_)

v = np.zeros((1,3))
v[0][2] = 1

T21 = np.dot(np.dot(K2, U_), v.T)
T22 = np.dot(np.dot(-K2, U_), v.T)

P21 = np.concatenate((R21, T21), axis=1)
P22 = np.concatenate((R21, T22), axis=1)
P23 = np.concatenate((R22, T21), axis=1)
P24 = np.concatenate((R22, T22), axis=1)

init_printing(wrap_line=False)

pprint(Matrix(P1))
print()
pprint(Matrix(P21))
print()
pprint(Matrix(P22))
print()
pprint(Matrix(P23))
print()
pprint(Matrix(P24))
