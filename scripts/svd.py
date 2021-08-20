import numpy as np
from sympy import *
import math

# A = [
# [2.267e+00,  6.956e-01,  1.550e+00,  9.485e-01,  2.910e-01,  6.485e-01,  1.463e+00,  4.487e-01,  1.000e+00],
# [2.856e+00,  -1.198e+00, -1.761e+00, -7.071e-01, 2.966e-01,  4.360e-01,  -1.622e+00, 6.804e-01,  1.000e+00],
# [5.674e-01,  -3.737e-01, -6.394e-01, -4.198e-01, 2.764e-01,  4.730e-01,  -8.874e-01, 5.844e-01,  1.000e+00],
# [8.901e-01,  4.129e-01,  8.505e-01,  6.378e-01,  2.959e-01,  6.094e-01,  1.047e+00,  4.855e-01,  1.000e+00],
# [2.267e+00,  -1.889e-01, 1.550e+00,  -4.624e-02, 3.852e-03,  -3.162e-02, 1.463e+00,  -1.218e-01, 1.000e+00],
# [2.856e+00,  1.219e+00,  -1.761e+00, 1.154e+00,  4.928e-01,  -7.117e-01, -1.622e+00, -6.924e-01, 1.000e+00],
# [5.674e-01,  4.427e-01,  -6.394e-01, 6.316e-01,  4.928e-01,  -7.117e-01, -8.874e-01, -6.924e-01, 1.000e+00],
# [8.901e-01,  -5.888e-01, 8.505e-01,  -7.449e-01, 4.928e-01,  -7.117e-01, 1.047e+00,  -6.924e-01, 1.000e+00],
# ]

# A = Matrix(A)

# print(A.nullspace())

A = [
 402037,  402037,     650,  371111,  371111,     600, 618.518, 618.518,       1,
 125301,  227820, 379.699,  204699,  372180, 620.301,     330,     600,       1,
 189941,  267866, 444.056,  261730,  369107, 611.888, 427.742, 603.226,       1,
 342086,  377398, 613.548,  336331,  371049, 603.226, 557.554, 615.108,       1,
 402037,  363518,     650,  340185,  307593,     550, 618.518, 559.259,       1,
 125301,  189850, 379.699,  165000,  250000,     500,     330,     500,       1,
 189941,  222028, 444.056,  213871,  250000,     500, 427.742,     500,       1,
 342086,  306774, 613.548,  278777,  250000,     500, 557.554,     500,       1,
]

A = np.reshape(A, (8, 9))
B = np.dot(A.T, A)
S,V,D = np.linalg.svd(B)

print(S)
print()
print(V)
print()
print(D)

v = D[:,8]

print(v)
print(np.dot(A, v))

exit()

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
