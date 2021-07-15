import math

n = 4

x = [0, 3, 1, 7, 0, 4, 1, 6, 3, 0]

out = [0, 0, 0, 0, 0, 0, 0, 0, 0]


# 0 in[0] -> tmp[0]
# 1 in[2] -> tmp[2]
# 2 in[4] -> tmp[4]
# 3 in[6] -> tmp[6]


# 0 in[1] -> tmp[1]
# 1 in[3] -> tmp[3]
# 2 in[5] -> tmp[5]
# 3 in[7] -> tmp[7]
# for lid in range(1):
#     stride = 1
#     while stride <= 4:
#         index = stride*(lid + 1) * 2 - 1;
#         if(index < 4):
#             print(index, " ", index - stride)
#         stride <<= 1
    


offset = 1
d = n>>1
while d > 0:
    for lid in range(3):
        if lid < d:
            ai = offset*(2*lid+1) - 1
            bi = offset*(2*lid+2) - 1
            print(ai, bi)
    offset *= 2
    d >>= 1


