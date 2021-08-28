import numpy as np
from numpy.lib.twodim_base import diag
import math
import cv2

pts1 = [
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
]

pts2 = [
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
]

pts1 = np.reshape(pts1, (8, 2))
pts2 = np.reshape(pts2, (8, 2))

E = [
  1.22081e-07,   5.12000e+00,   2.59271e-07, 
  5.12000e+00,  -3.81352e-07,   5.12000e+00, 
 -1.22081e-07,  -5.12000e+00,  -2.59270e-07,
]

E = np.reshape(E, (3,3))

def estimateRotationAndTranslation(E):
    U,D,V_T = np.linalg.svd(E, full_matrices=True)

    e = (D[0] + D[1])/2

    D[0] = e
    D[1] = e
    D[2] = 0

    E = np.dot(np.dot(U, np.diag(D)), V_T)

    U, _, V_T = np.linalg.svd(E)
    # U = -1*U
    # V_T = -1*V_T
    # print("V_T")
    # print(V_T)
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

    return (R1, R2, t, -t)

def triangulate(p, p_, P, P_):
    u = p[0] 
    u_ = p_[0] 
    v = p[1] 
    v_ = p_[1]
 
    A = [
        u*P[2][0] - P[0][0],  u*P[2][1] - P[0][1], u*P[2][2] - P[0][2], u*P[2][3] - P[0][3],
        v*P[2][0] - P[1][0],  v*P[2][1] - P[1][1], v*P[2][2] - P[1][2], v*P[2][3] - P[1][3],
        u_*P_[2][0] - P_[0][0],  u_*P_[2][1] - P_[0][1], u_*P_[2][2] - P_[0][2], u_*P_[2][3] - P_[0][3],
        v_*P_[2][0] - P_[1][0],  v_*P_[2][1] - P_[1][1], v_*P_[2][2] - P_[1][2], v_*P_[2][3] - P_[1][3],
    ]

    a = [
        u*P[2][0] - P[0][0],  u*P[2][1] - P[0][1], u*P[2][2] - P[0][2],
        v*P[2][0] - P[1][0],  v*P[2][1] - P[1][1], v*P[2][2] - P[1][2],
        u_*P_[2][0] - P_[0][0],  u_*P_[2][1] - P_[0][1], u_*P_[2][2] - P_[0][2],
        v_*P_[2][0] - P_[1][0],  v_*P_[2][1] - P_[1][1], v_*P_[2][2] - P_[1][2],
    ]

    a = np.reshape(a, (4, 3))

    b = [
        -(u*P[2][3] - P[0][3]),
        -(v*P[2][3] - P[1][3]),
        -(u_*P_[2][3] - P_[0][3]),
        -(v_*P_[2][3] - P_[1][3]),
    ]
    b = np.reshape(b, (4, 1))
    x = np.zeros((3,1))

    r = cv2.solve(a, b, x, cv2.DECOMP_SVD);
    print("(",r[1][0][0],",", r[1][1][0],",", r[1][2][0], ")")
    # A = np.reshape(A, (4,4))
    # S,D,V_T =  np.linalg.svd(A)
    # print("D")
    # print(D)
    # print(V_T)


def chooseRealizableSolution(Rotations, Translations, K0, K1, pts1, pts2):
    cam0 = np.dot(K0, np.concatenate((np.eye(3), np.zeros((1,3)).T), axis=1))

    cam1 = np.dot(K1, np.concatenate((Rotations[0].T, Translations[0].T), axis=1))
    for i in range(8):
        triangulate([pts1[i][0], pts1[i][1]], [pts2[i][0], pts2[i][1]], cam0, cam1)
    print()

    cam1 = np.dot(K1, np.concatenate((Rotations[0].T, Translations[1].T), axis=1))
    for i in range(8):
        triangulate([pts1[i][0], pts1[i][1]], [pts2[i][0], pts2[i][1]], cam0, cam1)
    print()

    cam1 = np.dot(K1, np.concatenate((Rotations[1].T, Translations[0].T), axis=1))
    for i in range(8):
        triangulate([pts1[i][0], pts1[i][1]], [pts2[i][0], pts2[i][1]], cam0, cam1)
    print()

    cam1 = np.dot(K1, np.concatenate((Rotations[1].T, Translations[1].T), axis=1))
    for i in range(8):
        triangulate([pts1[i][0], pts1[i][1]], [pts2[i][0], pts2[i][1]], cam0, cam1)
    print()



def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6



def rotationMatrixToEulerAngles(R) :
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])

    singular = sy < 1e-6

    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])
    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0

    print(np.round(np.rad2deg(x)), np.round(np.rad2deg(y)), np.round(np.rad2deg(z)))

K = [1600, 0, 500, 0, 1600, 500, 0, 0, 1]
K = np.reshape(K, (3,3))

R1, R2, t1, t2 = estimateRotationAndTranslation(E)

print(t1)
print(t2)
print(R1)
print(R2)

rotationMatrixToEulerAngles(R1)
rotationMatrixToEulerAngles(R2)
rotationMatrixToEulerAngles(R2)

R4 = [ 1.88248688e-08,  1.46553336e-05, -1.00000000e+00,
 1.21329510e-05,  1.00000000e+00,  1.46553338e-05,
 1.00000000e+00, -1.21329513e-05,  1.86470562e-08]
R4 = np.reshape(R4,(3,3))

print(np.linalg.det(R4))

R3 = [-7.85501760e-04,  8.49057726e-03, -9.99963646e-01,
  8.47538912e-03,  9.99928095e-01,  8.48361773e-03,
  9.99963775e-01, -8.46841711e-03, -8.57406225e-04]

R3 = np.reshape(R3, (3,3))
print(np.linalg.det(R3))
print("this")
rotationMatrixToEulerAngles(R4)
rotationMatrixToEulerAngles(R3)

# chooseRealizableSolution([R1,R2], [t1, t2], K, K, pts1, pts2)

# t = [ 
#     7.01373374e-01,
#     -2.95597654e-04,
#     7.12794011e-01
# ]
# t = np.reshape(t, (1,3))

# cam0 = np.dot(K, np.concatenate((np.eye(3), np.zeros((1,3)).T), axis=1))

# cam1 = np.dot(K, np.concatenate((R2.T, t.T), axis=1))

# for i in range(8):
#     triangulate([pts1[i][0], pts1[i][1]], [pts2[i][0], pts2[i][1]], cam0, cam1)
# print()

# cam1 = np.dot(K, np.concatenate((R2.T, -t.T), axis=1))
# for i in range(8):
#     triangulate([pts1[i][0], pts1[i][1]], [pts2[i][0], pts2[i][1]], cam0, cam1)
# print()

# cam1 = np.dot(K, np.concatenate((R2, -t.T), axis=1))
# for i in range(8):
#     triangulate([pts1[i][0], pts1[i][1]], [pts2[i][0], pts2[i][1]], cam0, cam1)
# print()
