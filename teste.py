import math

def map():
    i = 3
    j = 4
    k = 5
    
    t = i * j * k
    for p in range(0, t):

        print('%.2f'%(math.floor((math.floor((p / k))/j))%i), '%.2f'%(math.floor((p / k))%j), '%.2f'%(p % k))

map()

