from sympy import *

h = Matrix([Symbol("h1"), Symbol("h2"), Symbol("h3"), Symbol("h4"),  Symbol("h5"),  Symbol("h6"),  Symbol("h7"),  Symbol("h8"),  Symbol("h9")])

B = Matrix([
    [Symbol("p1[0]"), Symbol("p1[1]"), 1, 0,0,0,-1*Symbol("p1[0]")*Symbol("p2[0]"), -1*Symbol("p1[1]")*Symbol("p2[0]"), 0],
    [0,0,0,Symbol("p1[0]"), Symbol("p1[1]"), 1, -1*Symbol("p1[0]")*Symbol("p2[1]"), -1*Symbol("p1[1]")*Symbol("p2[1]"), 0],

    [Symbol("p3[0]"), Symbol("p3[1]"), 1, 0,0,0,-1*Symbol("p3[0]")*Symbol("p4[0]"), -1*Symbol("p3[1]")*Symbol("p4[0]"), 0],
    [0,0,0,Symbol("p3[0]"), Symbol("p3[1]"), 1, -1*Symbol("p3[0]")*Symbol("p4[1]"), -1*Symbol("p3[1]")*Symbol("p4[1]"), 0],

    [Symbol("p5[0]"), Symbol("p5[1]"), 1, 0,0,0,-1*Symbol("p5[0]")*Symbol("p6[0]"), -1*Symbol("p5[1]")*Symbol("p6[0]"), 0],
    [0,0,0,Symbol("p5[0]"), Symbol("p5[1]"), 1, -1*Symbol("p5[0]")*Symbol("p6[1]"), -1*Symbol("p5[1]")*Symbol("p6[1]"), 0],
])

nullspace = B.nullspace()

A = nullspace[0]
B = nullspace[1]
C = nullspace[2]

# p11[0] = u1,1
# p11[1] = v1,1
# p21[0] = u2,1
# p21[1] = v2,1

# p12[0] = u1,2
# p12[1] = v1,2
# p22[0] = u2,2
# p22[1] = v2,2

beta0 = -cos(Symbol("a[1]"))*Symbol("c[0]")*Symbol("d[6]")*Symbol("p22[1]")*sin(Symbol("a[0]"))
beta0 += cos(Symbol("a[1]"))*Symbol("c[3]")*Symbol("d[6]")*Symbol("p21[0]")*sin(Symbol("a[0]"))
beta0 += cos(Symbol("a[1]"))*Symbol("c[6]")*Symbol("d[0]")*Symbol("p22[1]")*sin(Symbol("a[0]"))

beta0 -= cos(Symbol("a[1]"))*Symbol("c[6]")*Symbol("d[3]")*Symbol("p21[0]")*sin(Symbol("a[0]"))
beta0 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("c[3]")*Symbol("d[6]")*Symbol("p21[1]")
beta0 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("c[3]")*Symbol("d[6]")*Symbol("p22[1]")

beta0 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("c[6]")*Symbol("d[3]")*Symbol("p21[1]")
beta0 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("c[6]")*Symbol("d[3]")*Symbol("p22[1]")
beta0 -= Symbol("c[0]")*Symbol("d[6]")*Symbol("p21[0]")*sin(Symbol("a[1]"))*sin(Symbol("a[0]"))

beta0 += Symbol("c[0]")*Symbol("d[6]")*Symbol("p22[0]")*sin(Symbol("a[1]"))*sin(Symbol("a[0]"))
beta0 += Symbol("c[6]")*Symbol("d[0]")*Symbol("p21[0]")*sin(Symbol("a[1]"))*sin(Symbol("a[0]"))
beta0 -= Symbol("c[6]")*Symbol("d[0]")*Symbol("p22[0]")*sin(Symbol("a[1]"))*sin(Symbol("a[0]"))

beta0 += cos(Symbol("a[0]"))*Symbol("c[0]")*Symbol("d[6]")*Symbol("p21[1]")*sin(Symbol("a[1]"))
beta0 -= cos(Symbol("a[0]"))*Symbol("c[3]")*Symbol("d[6]")*Symbol("p22[0]")*sin(Symbol("a[1]"))
beta0 -= cos(Symbol("a[0]"))*Symbol("c[6]")*Symbol("d[0]")*Symbol("p21[1]")*sin(Symbol("a[1]"))

beta0 += cos(Symbol("a[0]"))*Symbol("c[6]")*Symbol("d[3]")*Symbol("p22[0]")*sin(Symbol("a[1]"))
beta0 += cos(Symbol("a[1]"))*Symbol("c[0]")*Symbol("d[3]")*sin(Symbol("a[0]"))
beta0 -= cos(Symbol("a[1]"))*Symbol("c[3]")*Symbol("d[0]")*sin(Symbol("a[0]"))

beta0 -= cos(Symbol("a[0]"))*Symbol("c[0]")*Symbol("d[3]")*sin(Symbol("a[1]"))
beta0 += cos(Symbol("a[0]"))*Symbol("c[3]")*Symbol("d[0]")*sin(Symbol("a[1]"))

beta1 = +cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("c[6]")*Symbol("p22[1]")*sin(Symbol("a[0]"))
beta1 += cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("p21[0]")*sin(Symbol("a[0]"))
beta1 += cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("c[0]")*Symbol("p22[1]")*sin(Symbol("a[0]"))

beta1 -= cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("p21[0]")*sin(Symbol("a[0]"))
beta1 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("p21[1]")
beta1 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("p22[1]")

beta1 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("p21[1]")
beta1 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("p22[1]")
beta1 -= Symbol("b[0]")*Symbol("c[6]")*Symbol("p21[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

beta1 += Symbol("b[0]")*Symbol("c[6]")*Symbol("p22[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
beta1 += Symbol("b[6]")*Symbol("c[0]")*Symbol("p21[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
beta1 -= Symbol("b[6]")*Symbol("c[0]")*Symbol("p22[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

beta1 += cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("c[6]")*Symbol("p21[1]")*sin(Symbol("a[1]"))
beta1 -= cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("p22[0]")*sin(Symbol("a[1]"))
beta1 -= cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[0]")*Symbol("p21[1]")*sin(Symbol("a[1]"))

beta1 += cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("p22[0]")*sin(Symbol("a[1]"))
beta1 += cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("c[3]")*sin(Symbol("a[0]"))
beta1 -= cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("c[0]")*sin(Symbol("a[0]"))

beta1 -= cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("c[3]")*sin(Symbol("a[1]"))
beta1 += cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[0]")*sin(Symbol("a[1]"))

beta = beta0/beta1

gama0 = -cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("d[6]")*Symbol("p22[1]")*sin(Symbol("a[0]"))
gama0 += cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("d[6]")*Symbol("p21[0]")*sin(Symbol("a[0]"))
gama0 += cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("d[0]")*Symbol("p22[1]")*sin(Symbol("a[0]"))

gama0 -= cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("d[3]")*Symbol("p21[0]")*sin(Symbol("a[0]"))
gama0 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("d[6]")*Symbol("p21[1]")
gama0 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("d[6]")*Symbol("p22[1]")

gama0 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("d[3]")*Symbol("p21[1]")
gama0 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("d[3]")*Symbol("p22[1]")
gama0 -= Symbol("b[0]")*Symbol("d[6]")*Symbol("p21[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

gama0 += Symbol("b[0]")*Symbol("d[6]")*Symbol("p22[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
gama0 += Symbol("b[6]")*Symbol("d[0]")*Symbol("p21[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
gama0 -= Symbol("b[6]")*Symbol("d[0]")*Symbol("p22[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

gama0 += cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("d[6]")*Symbol("p21[1]")*sin(Symbol("a[1]"))
gama0 -= cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("d[6]")*Symbol("p22[0]")*sin(Symbol("a[1]"))
gama0 -= cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("d[0]")*Symbol("p21[1]")*sin(Symbol("a[1]"))

gama0 += cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("d[3]")*Symbol("p22[0]")*sin(Symbol("a[1]"))
gama0 += cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("d[3]")*sin(Symbol("a[0]"))
gama0 -= cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("d[0]")*sin(Symbol("a[0]"))

gama0 -= cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("d[3]")*sin(Symbol("a[1]"))
gama0 += cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("d[0]")*sin(Symbol("a[1]"))

gama1 = -cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("c[6]")*Symbol("p22[1]")*sin(Symbol("a[0]"))
gama1 += cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("p21[0]")*sin(Symbol("a[0]"))
gama1 += cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("c[0]")*Symbol("p22[1]")*sin(Symbol("a[0]"))

gama1 -= cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("p21[0]")*sin(Symbol("a[0]"))
gama1 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("p21[1]")
gama1 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("p22[1]")

gama1 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("p21[1]")
gama1 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("p22[1]")
gama1 -= Symbol("b[0]")*Symbol("c[6]")*Symbol("p21[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

gama1 += Symbol("b[0]")*Symbol("c[6]")*Symbol("p22[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
gama1 += Symbol("b[6]")*Symbol("c[0]")*Symbol("p21[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
gama1 -= Symbol("b[6]")*Symbol("c[0]")*Symbol("p22[0]")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

gama1 += cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("c[6]")*Symbol("p21[1]")*sin(Symbol("a[1]"))
gama1 -= cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("p22[0]")*sin(Symbol("a[1]"))
gama1 -= cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[0]")*Symbol("p21[1]")*sin(Symbol("a[1]"))

gama1 += cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("p22[0]")*sin(Symbol("a[1]"))
gama1 += cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("c[3]")*sin(Symbol("a[0]"))
gama1 -= cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("c[0]")*sin(Symbol("a[0]"))

gama1 -= cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("c[3]")*sin(Symbol("a[1]"))
gama1 += cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[0]")*sin(Symbol("a[1]"))

gama = -gama0 / gama1



# p11[0] = u1,1
# p11[1] = v1,1
# p21[0] = u2,1
# p21[1] = v2,1

# p12[0] = u1,2
# p12[1] = v1,2
# p22[0] = u2,2
# p22[1] = v2,2

# D = Matrix([
#     [Symbol("u1,1")*Symbol("u2,1"), Symbol("v1,1")*Symbol("u2,1"), Symbol("u2,1"), Symbol("u1,1")*Symbol("v2,1"), Symbol("v1,1")*Symbol("v  2,1"), Symbol("v2,1"), Symbol("u1,1"), Symbol("v1,1"), 1],
#     [Symbol("u1,2")*Symbol("u2,2"), Symbol("v1,2")*Symbol("u2,2"), Symbol("u2,2"), Symbol("u1,2")*Symbol("v2,2"), Symbol("v1,2")*Symbol("v2,2"), Symbol("v2,2"), Symbol("u1,2"), Symbol("v1,2"), 1],
#     [Symbol("u1,3")*Symbol("u2,3"), Symbol("v1,3")*Symbol("u2,3"), Symbol("u2,3"), Symbol("u1,3")*Symbol("v2,3"), Symbol("v1,3")*Symbol("v2,3"), Symbol("v2,3"), Symbol("u1,3"), Symbol("v1,3"), 1],
#     [Symbol("u1,4")*Symbol("u2,4"), Symbol("v1,4")*Symbol("u2,4"), Symbol("u2,4"), Symbol("u1,4")*Symbol("v2,4"), Symbol("v1,4")*Symbol("v2,4"), Symbol("v2,4"), Symbol("u1,4"), Symbol("v1,4"), 1],
#     [Symbol("u1,5")*Symbol("u2,5"), Symbol("v1,5")*Symbol("u2,5"), Symbol("u2,5"), Symbol("u1,5")*Symbol("v2,5"), Symbol("v1,5")*Symbol("v2,5"), Symbol("v2,5"), Symbol("u1,5"), Symbol("v1,5"), 1],
#     [Symbol("u1,6")*Symbol("u2,6"), Symbol("v1,6")*Symbol("u2,6"), Symbol("u2,6"), Symbol("u1,6")*Symbol("v2,6"), Symbol("v1,6")*Symbol("v2,6"), Symbol("v2,6"), Symbol("u1,6"), Symbol("v1,6"), 1],
#     [Symbol("u1,7")*Symbol("u2,7"), Symbol("v1,7")*Symbol("u2,7"), Symbol("u2,7"), Symbol("u1,7")*Symbol("v2,7"), Symbol("v1,7")*Symbol("v2,7"), Symbol("v2,7"), Symbol("u1,7"), Symbol("v1,7"), 1],
# ])

# D_nullspace = D.nullspace()

# print(D_nullspace)

x = Symbol("x")
n = 1-x

e = Matrix([Symbol("e1"), Symbol("e2"), Symbol("e3"), Symbol("e4"), Symbol("e5"), Symbol("e6"), Symbol("e7"), Symbol("e8"), Symbol("e9")])
g = Matrix([Symbol("g1"), Symbol("g2"), Symbol("g3"), Symbol("g4"), Symbol("g5"), Symbol("g6"), Symbol("g7"), Symbol("g8"), Symbol("g9")])

f = Matrix(x*e + n*g)

F = Matrix([
    [f[0], f[1], f[2]],
    [f[3], f[4], f[5]],
    [f[6], f[7], f[8]],
])

d = det(F)
d = Poly(expand(d), Symbol("x"))
coeffs = d.coeffs()
print(factor(d))
print(coeffs[0])
print(coeffs[1])
print(coeffs[2])
print(coeffs[3])

# print(len(d.coeffs()))
# rts = roots(d)
# print(rts)
f = open("src/bundle/codegen/Homography.cpp", "a")

f.truncate(0)

f.write("#include <cmath>\n")
f.write("#include \"Homography.hpp\"\n")


f.write("double beta(double* b, double* c, double* d, double* p11, double* p21, double* p12, double* p22, double* a)\n{\n return %s;\n}\n\n" % cxxcode(beta))
f.write("double gamma(double* b, double* c, double* d, double* p11, double* p21, double* p12, double* p22, double* a)\n{\n return %s;\n}\n\n" % cxxcode(gama))


f.write("double nullSpaceVecB0(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(A[0]))
f.write("double nullSpaceVecB1(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(A[1]))
f.write("double nullSpaceVecB2(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(A[2]))
f.write("double nullSpaceVecB3(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(A[3]))
f.write("double nullSpaceVecB4(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(A[4]))
f.write("double nullSpaceVecB5(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(A[5]))
f.write("double nullSpaceVecB6(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(A[6]))
f.write("double nullSpaceVecB7(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(A[7]))
f.write("double nullSpaceVecB8(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(A[8]))

f.write("double nullSpaceVecC0(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(B[0]))
f.write("double nullSpaceVecC1(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(B[1]))
f.write("double nullSpaceVecC2(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(B[2]))
f.write("double nullSpaceVecC3(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(B[3]))
f.write("double nullSpaceVecC4(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(B[4]))
f.write("double nullSpaceVecC5(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(B[5]))
f.write("double nullSpaceVecC6(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(B[6]))
f.write("double nullSpaceVecC7(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(B[7]))
f.write("double nullSpaceVecC8(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(B[8]))

f.write("double nullSpaceVecD0(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(C[0]))
f.write("double nullSpaceVecD1(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(C[1]))
f.write("double nullSpaceVecD2(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(C[2]))
f.write("double nullSpaceVecD3(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(C[3]))
f.write("double nullSpaceVecD4(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(C[4]))
f.write("double nullSpaceVecD5(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(C[5]))
f.write("double nullSpaceVecD6(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(C[6]))
f.write("double nullSpaceVecD7(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(C[7]))
f.write("double nullSpaceVecD8(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)\n{\n return %s;\n}\n\n" % cxxcode(C[8]))

f.close()


f = open("src/bundle/codegen/Homography.hpp", "a")

f.truncate(0)

f.write("#pragma once\n")

f.write("double beta(double* b, double* c, double* d, double* p11, double* p21, double* p12, double* p22, double* a);\n")
f.write("double gamma(double* b, double* c, double* d, double* p11, double* p21, double* p12, double* p22, double* a);\n")

f.write("double nullSpaceVecB0(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecB1(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecB2(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecB3(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecB4(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecB5(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecB6(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecB7(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecB8(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")

f.write("double nullSpaceVecC0(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecC1(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecC2(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecC3(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecC4(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecC5(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecC6(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecC7(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecC8(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")

f.write("double nullSpaceVecD0(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecD1(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecD2(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecD3(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecD4(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecD5(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecD6(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecD7(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")
f.write("double nullSpaceVecD8(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6);\n")

f.close()


# print(nullspace[0][0])

