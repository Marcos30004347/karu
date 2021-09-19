import math

def map():
    i = 3
    j = 4
    k = 5
    
    t = i * j * k
    for idx in range(0, t):
        print(
            '%.2f'%(math.floor((math.floor((idx / k))/j))%i),
            '%.2f'%(math.floor((idx / k))%j),
            '%.2f'%(idx % k)
        )

map()

