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


#  [u[0],[image], v[0],[image], 1], [u[1],[image], v[1],[image]]
beta0 = -1*cos(Symbol("a[1]"))*Symbol("c[0]")*Symbol("d[6]")*Symbol("v2,2")*sin(Symbol("a[0]"))
beta0 += cos(Symbol("a[1]"))*Symbol("c[3]")*Symbol("d[6]")*Symbol("u2,1")*sin(Symbol("a[0]"))
beta0 += cos(Symbol("a[1]"))*Symbol("c[6]")*Symbol("d[0]")*Symbol("v2,2")*sin(Symbol("a[0]"))

beta0 -= cos(Symbol("a[1]"))*Symbol("c[6]")*Symbol("d[3]")*Symbol("u2,1")*sin(Symbol("a[0]"))
beta0 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("c[3]")*Symbol("d[6]")*Symbol("v2,1")
beta0 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("c[3]")*Symbol("d[6]")*Symbol("v2,2")

beta0 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("c[6]")*Symbol("d[3]")*Symbol("v2,1")
beta0 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("c[6]")*Symbol("d[3]")*Symbol("v2,2")
beta0 -= Symbol("c[0]")*Symbol("d[6]")*Symbol("u2,1")*sin(Symbol("a[1]"))*sin(Symbol("a[0]"))

beta0 += Symbol("c[0]")*Symbol("d[6]")*Symbol("u2,2")*sin(Symbol("a[1]"))*sin(Symbol("a[0]"))
beta0 += Symbol("c[6]")*Symbol("d[0]")*Symbol("u2,1")*sin(Symbol("a[1]"))*sin(Symbol("a[0]"))
beta0 -= Symbol("c[6]")*Symbol("d[0]")*Symbol("u2,2")*sin(Symbol("a[1]"))*sin(Symbol("a[0]"))

beta0 += cos(Symbol("a[0]"))*Symbol("c[0]")*Symbol("d[6]")*Symbol("v2,1")*sin(Symbol("a[1]"))
beta0 -= cos(Symbol("a[0]"))*Symbol("c[3]")*Symbol("d[6]")*Symbol("u2,2")*sin(Symbol("a[1]"))
beta0 -= cos(Symbol("a[0]"))*Symbol("c[6]")*Symbol("d[0]")*Symbol("v2,1")*sin(Symbol("a[1]"))

beta0 += cos(Symbol("a[0]"))*Symbol("c[6]")*Symbol("d[3]")*Symbol("u2,2")*sin(Symbol("a[1]"))
beta0 += cos(Symbol("a[1]"))*Symbol("c[0]")*Symbol("d[3]")*sin(Symbol("a[0]"))
beta0 -= cos(Symbol("a[1]"))*Symbol("c[3]")*Symbol("d[0]")*sin(Symbol("a[0]"))

beta0 -= cos(Symbol("a[0]"))*Symbol("c[0]")*Symbol("d[3]")*sin(Symbol("a[1]"))
beta0 += cos(Symbol("a[0]"))*Symbol("c[3]")*Symbol("d[0]")*sin(Symbol("a[1]"))

beta1 = cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("c[6]")*Symbol("v2,2")*sin(Symbol("a[0]"))
beta1 += cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("u2,1")*sin(Symbol("a[0]"))
beta1 += cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("c[0]")*Symbol("v2,2")*sin(Symbol("a[0]"))

beta1 -= cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("u2,1")*sin(Symbol("a[0]"))
beta1 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("v2,1")
beta1 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("v2,2")

beta1 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("v2,1")
beta1 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("v2,2")
beta1 -= Symbol("b[0]")*Symbol("c[6]")*Symbol("u2,1")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

beta1 += Symbol("b[0]")*Symbol("c[6]")*Symbol("u2,2")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
beta1 += Symbol("b[6]")*Symbol("c[0]")*Symbol("u2,1")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
beta1 -= Symbol("b[6]")*Symbol("c[0]")*Symbol("u2,2")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

beta1 += cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("c[6]")*Symbol("v2,1")*sin(Symbol("a[1]"))
beta1 -= cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("u2,2")*sin(Symbol("a[1]"))
beta1 -= cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[0]")*Symbol("v2,1")*sin(Symbol("a[1]"))

beta1 += cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("u2,2")*sin(Symbol("a[1]"))
beta1 += cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("c[3]")*sin(Symbol("a[0]"))
beta1 -= cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("c[0]")*sin(Symbol("a[0]"))

beta1 -= cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("c[3]")*sin(Symbol("a[1]"))
beta1 -= cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[0]")*sin(Symbol("a[1]"))

beta = beta0/beta1



gama0 = -1*cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("d[6]")*Symbol("v2,2")*sin(Symbol("a[0]"))
gama0 += cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("d[6]")*Symbol("u2,1")*sin(Symbol("a[0]"))
gama0 += cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("d[0]")*Symbol("v2,2")*sin(Symbol("a[0]"))

gama0 -= cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("d[3]")*Symbol("u2,1")*sin(Symbol("a[0]"))
gama0 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("d[6]")*Symbol("v2,1")
gama0 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("d[6]")*Symbol("v2,2")

gama0 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("d[3]")*Symbol("v2,1")
gama0 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("d[3]")*Symbol("v2,2")
gama0 -= Symbol("b[0]")*Symbol("d[6]")*Symbol("u2,1")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

gama0 += Symbol("b[0]")*Symbol("d[6]")*Symbol("u2,2")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
gama0 += Symbol("b[6]")*Symbol("d[0]")*Symbol("u2,1")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
gama0 -= Symbol("b[6]")*Symbol("d[0]")*Symbol("u2,2")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

gama0 += cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("d[6]")*Symbol("v2,1")*sin(Symbol("a[1]"))
gama0 -= cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("d[6]")*Symbol("v2,2")*sin(Symbol("a[1]"))
gama0 += cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("d[0]")*Symbol("v2,1")*sin(Symbol("a[1]"))

gama0 += cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("d[3]")*Symbol("u2,2")*sin(Symbol("a[1]"))
gama0 += cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("d[3]")*sin(Symbol("a[0]"))
gama0 += cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("d[0]")*sin(Symbol("a[0]"))

gama0 -= cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("d[3]")*sin(Symbol("a[1]"))
gama0 += cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("d[0]")*sin(Symbol("a[1]"))

gama1 = -1*cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("c[6]")*Symbol("v2,2")*sin(Symbol("a[0]"))
gama1 += cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("u2,1")*sin(Symbol("a[0]"))
gama1 += cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("c[0]")*Symbol("v2,2")*sin(Symbol("a[0]"))

gama1 -= cos(Symbol("a[1]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("u2,1")*sin(Symbol("a[0]"))
gama1 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("v2,1")
gama1 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("v2,2")

gama1 += cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("v2,1")
gama1 -= cos(Symbol("a[1]"))*cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("v2,2")
gama1 -= Symbol("b[0]")*Symbol("c[6]")*Symbol("u2,1")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

gama1 += Symbol("b[0]")*Symbol("c[6]")*Symbol("u2,2")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
gama1 += Symbol("b[6]")*Symbol("c[0]")*Symbol("u2,1")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))
gama1 -= Symbol("b[6]")*Symbol("c[0]")*Symbol("u2,2")*sin(Symbol("a[0]"))*sin(Symbol("a[1]"))

gama1 += cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("c[6]")*Symbol("v2,1")*sin(Symbol("a[1]"))
gama1 -= cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[6]")*Symbol("u2,2")*sin(Symbol("a[1]"))
gama1 -= cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[0]")*Symbol("v2,1")*sin(Symbol("a[1]"))

gama1 += cos(Symbol("a[0]"))*Symbol("b[6]")*Symbol("c[3]")*Symbol("u2,2")*sin(Symbol("a[1]"))
gama1 += cos(Symbol("a[1]"))*Symbol("b[0]")*Symbol("c[3]")*sin(Symbol("a[0]"))
gama1 -= cos(Symbol("a[1]"))*Symbol("b[3]")*Symbol("c[0]")*sin(Symbol("a[0]"))

gama1 -= cos(Symbol("a[0]"))*Symbol("b[0]")*Symbol("c[3]")*sin(Symbol("a[1]"))
gama1 += cos(Symbol("a[0]"))*Symbol("b[3]")*Symbol("c[0]")*sin(Symbol("a[1]"))

gama = -1*gama0 / gama1




f = open("src/bundle/codegen/Homography.cpp", "a")

f.truncate(0)

f.write("#include <cmath>\n")


f.write("double beta(double* c, double* b, double* d, double* p1, double* p2, double* p3, double* p4, double* a)\n{\n return %s;\n}\n\n" % cxxcode(beta))
f.write("double gamma(double* c, double* b, double* d, double* p1, double* p2, double* p3, double* p4, double* a)\n{\n return %s;\n}\n\n" % cxxcode(gama))


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

f.write("double beta(double* c, double* b, double* d, double* p1, double* p2, double* p3, double* p4, double* a);\n")
f.write("double gamma(double* c, double* b, double* d, double* p1, double* p2, double* p3, double* p4, double* a);\n")

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

