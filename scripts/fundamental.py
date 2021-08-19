import numpy as np
import cv2

pts1 = [
7.407e-02,
7.407e-02,
-1.062e-01,
6.250e-02,
-4.516e-02,
6.452e-02,
3.597e-02,
7.194e-02,
7.407e-02,
-2.679e-08,
-1.063e-01,
-2.679e-08,
-4.516e-02,
-2.679e-08,
3.597e-02,
-2.679e-08,
]
pts1 = np.reshape(pts1,(8,2))

pts2 = [
9.375e-02,
6.250e-02,
-7.519e-02,
7.519e-02,
-3.497e-02,
6.993e-02,
7.097e-02,
6.452e-02,
9.375e-02,
2.161e-08,
-7.519e-02,
2.564e-08,
-3.497e-02,
2.468e-08,
7.097e-02,
2.215e-08,
]

pts2 = np.reshape(pts2,(8,2))

K = [1, 0, 0, 0, 1, 0, 0, 0, 1]
K = np.reshape(K, (3,3))

print(K)

F, mask = cv2.findFundamentalMat(pts1,pts2,cv2.FM_8POINT, 3, 0.99)
# E, mask = cv2.findEssentialMat(pts1,pts2, K,cv2.FM_8POINT, 0.99)




E = np.dot(np.dot(K.T, F), K)

pts, r, t, mask = cv2.recoverPose(E,pts1, pts2, K)
r = np.linalg.inv(r)
t = t * -1

print(r)
print(t)
print(pts)
