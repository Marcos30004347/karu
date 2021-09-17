-- fill x with random numbers --
algorithm fill random 'x':
    if type(x) == 'matrix':
        for i in 1 to x.rows, j in 1 to x.cols: 
            x[i,j] = random()
        end
    end
    
    if type(x) == number:
        x = random()
    end

    return x
end

-- activation function --
a(v, t, k) = 3 * v

-- network function --
n(v, A, B, C, t, k) = a(a(a(v * A, t, k) * B, t, k) * C, t, k)

-- loss function --
l(outp, data) = (data - outp)^2

-- network parameters --
A = fill random 'matrix:3,3'
B = fill random 'matrix:3,3'
C = fill random 'matrix:3,3'

t = random()
k = random()

-- train network --
algorithm train network with input 'input' and desired output 'output':
    alpha := 0.1

    for epoch in 1 to 100, i in 1 to len(data):
        
        out = n(input[i], A, B, C)
        err = l(out, output[i])

        A -= alpha * diff(err, A)
        B -= alpha * diff(err, B)
        C -= alpha * diff(err, C)
        t -= alpha * diff(err, t)
        k -= alpha * diff(err, k)
    end
end

inp_test_data = [[1,2,3], [3,4,5], [6,7,8]]
out_test_data = [[1,2,3], [3,4,5], [6,7,8],]

train network with input 'inp_test_data' and desired output 'out_test_data'


-- some ideias aboud indexing --

a = tensor:10,10,10

a:(1 to 4),1,1 = 4

a:i,j,k = 4

a:1 = 2

a:2 =4

x = 4

A = matrix:3,3 [[1,2,3], [3,4,5], [6,7,8]]
B = matrix:3,3 [[1,2,3], [3,4,5], [6,7,8]]

B:1,1 := 3
B:2,1 := 3

D = tensor:3,3,1[[[1],[2],[3]], [[3],[4],[5]], [[6],[7],[8]]]

i = 3
j = 3
k = 1

D := tensor:i,j,k[[[1],[2],[3]], [[3],[4],[5]], [[6],[7],[8]]]


