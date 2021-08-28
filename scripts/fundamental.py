import numpy as np
import cv2

pts1 = np.array([
    6.185e+02,
    6.185e+02,
    3.300e+02,
    6.000e+02,
    4.277e+02,
    6.032e+02,
    5.576e+02,
    6.151e+02,
    6.185e+02,
    5.593e+02,
    3.300e+02,
    5.000e+02,
    4.277e+02,
    5.000e+02,
    5.576e+02,
    5.000e+02,
])

pts2 = np.array([
    6.500e+02,
    6.000e+02,
    3.797e+02,
    6.203e+02,
    4.441e+02,
    6.119e+02,
    6.135e+02,
    6.032e+02,
    6.500e+02,
    5.500e+02,
    3.797e+02,
    5.000e+02,
    4.441e+02,
    5.000e+02,
    6.135e+02,
    5.000e+02,
])

pts1 = np.reshape(pts1, (8, 2))
pts2 = np.reshape(pts2, (8, 2))

K = [1600, 0, 500, 0, 1600, 500, 0, 0, 1]
K = np.reshape(K, (3,3))

F, mask = cv2.findFundamentalMat(pts1, pts2, cv2.FM_8POINT, 3, 0.99)

print("Fundamental")
print(F)

p1 = np.array([pts1[0][0], pts1[0][1], 1], dtype=object)
p2 = np.array([pts2[0][0], pts2[0][1], 1], dtype=object)
p1 = np.reshape(p1, (3,1))
p2 = np.reshape(p2, (3,1))

print("******************")
print(p1)
print(p2)
print(p1.T @ F @ p2)
print("******************")

E = np.dot(np.dot(K.T, F), K)
print("Essential")
print(E)
print("##########################")
pts, r, t, mask = cv2.recoverPose(E,pts1, pts2, K)

r = r
t = t

print(r)
print(t)
print(pts)
np.set_printoptions(linewidth=1200)


E = [
8.226e-05,       5.121e+00,       1.196e-07,
5.119e+00,       -2.389e-04,      5.120e+00,
1.279e-05,       -5.120e+00,      -1.237e-06,
]
E = np.reshape(E, (3,3))

pts, r, t, mask = cv2.recoverPose(E,pts1, pts2, K)

r = r
t = t
print("ASDASSAD")
print(r)
print(t)
print(pts)


uv = [
	[pts1[0][0], pts1[0][1], pts2[0][0], pts2[0][1]],
	[pts1[1][0], pts1[1][1], pts2[1][0], pts2[1][1]],
	[pts1[2][0], pts1[2][1], pts2[2][0], pts2[2][1]],
	[pts1[3][0], pts1[3][1], pts2[3][0], pts2[3][1]],
	[pts1[4][0], pts1[4][1], pts2[4][0], pts2[4][1]],
	[pts1[5][0], pts1[5][1], pts2[5][0], pts2[5][1]],
	[pts1[6][0], pts1[6][1], pts2[6][0], pts2[6][1]],
	[pts1[7][0], pts1[7][1], pts2[7][0], pts2[7][1]],
]

def calc_F(uvMat):
    A = np.zeros((len(uvMat),9))
    # img1 x' y' x y im2
    for i in range(len(uvMat)):
        A[i][0] = uvMat[i][0]*uvMat[i][2]
        A[i][1] = uvMat[i][1]*uvMat[i][2]
        A[i][2] = uvMat[i][2]
        A[i][3] = uvMat[i][0]*uvMat[i][3]
        A[i][4] = uvMat[i][1]*uvMat[i][3]
        A[i][5] = uvMat[i][3] 
        A[i][6] = uvMat[i][0]
        A[i][7] = uvMat[i][1]
        A[i][8] = 1.0  

    print(A)
    # print(A)
    u,s,v = np.linalg.svd(A)
    print(s)
    print(v)
    # v = -1*v

    # print(s)
    f_vec = v.transpose()[:,8]
    # print("f_vec = ", f_vec)
    f_hat = np.reshape(f_vec, (3,3))
    # print("Fmat = ", f_hat)

    # Enforce rank(F) = 2 
    s,v,d = np.linalg.svd(f_hat)
    f_hat = s @ np.diag([*v[:2], 0]) @ d

    return f_hat
F = calc_F(uv)

print("asdasdasdasdasdasdasdasdsd")
print(F/F[2,2])
print("asdasdasdasdasdasdasdasdsd")

E = np.dot(np.dot(K.T, calc_F(uv)), K)

print(E)

pts, r, t, mask = cv2.recoverPose(E,pts1, pts2, K)

r = r
t = t

print(r)
print(t)
print(pts)
