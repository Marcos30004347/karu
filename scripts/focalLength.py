from sympy import *

init_printing(wrap_line=False)

U11 = Symbol("U[0][0]")
U12 = Symbol("U[0][1]")
U13 = Symbol("U[0][2]")
U21 = Symbol("U[1][0]")
U22 = Symbol("U[1][1]")
U23 = Symbol("U[1][2]")
U31 = Symbol("U[2][0]")
U32 = Symbol("U[2][1]")
U33 = Symbol("U[2][2]")

V11 = Symbol("V[0][0]")
V12 = Symbol("V[0][1]")
V13 = Symbol("V[0][2]")
V21 = Symbol("V[1][0]")
V22 = Symbol("V[1][1]")
V23 = Symbol("V[1][2]")
V31 = Symbol("V[2][0]")
V32 = Symbol("V[2][1]")
V33 = Symbol("V[2][2]")
r = Symbol("r")
s = Symbol("s")
x = Symbol("x")



M1 = Matrix([
    [U11 * V13, U12 * V13, U13 * V13, r * U11 * V11  + s * U12 * V12],
    [U11 * V23, U12 * V23, U13 * V23, r * U11 * V21  + s * U12 * V22],
    [U21 * V13, U22 * V13, U23 * V13, r * U21 * V11  + s * U22 * V12],
    [U21 * V23, U22 * V23, U23 * V23, r * U21 * V21  + s * U22 * V22],
])

Mx = Matrix([
    [-s * U13 * V11, -r * U13*V12, r*U12*V12 + s*U11*V11, r*s*U13*V13],
    [-s * U13 * V21, -r * U13*V22, r*U12*V22 + s*U11*V21, r*s*U13*V23],
    [-s * U23 * V11, -r * U23*V12, r*U22*V12 + s*U21*V11, r*s*U23*V13],
    [-s * U23 * V21, -r * U23*V22, r*U22*V22 + s*U21*V21, r*s*U23*V23],
])
print("det")
pprint(det(M1))
print()
print("det")
pprint(det(Mx))
print()
print()
print()
print()

d = det(M1 - x*Mx)

poly = Poly(d, x)

# pprint(poly)
# print(poly.all_coeffs())
# print("coeff 0")
# col_exp = collect(poly, x)
# pprint(col_exp)
# print("coeff 0")
# pprint(poly.all_coeffs()[0])
# print("coeff 1")
# pprint(poly.all_coeffs()[1])
# print("coeff 2")
# pprint(poly.all_coeffs()[2])
# print("coeff 3")
# pprint(poly.all_coeffs()[3])
# print("coeff 4")
# print(poly.all_coeffs())
f = open("src/bundle/codegen/Estimation.hpp", "a")

f.truncate(0)

f.write("#include <cmath>\n")

f.write("double polyCoeff3(double** U, double** V, double s, double r);\n")
f.write("double polyCoeff2(double** U, double** V, double s, double r);\n")
f.write("double polyCoeff1(double** U, double** V, double s, double r);\n")
f.write("double polyCoeff0(double** U, double** V, double s, double r);\n")

f.close()

f = open("src/bundle/codegen/Estimation.cpp", "a")

f.truncate(0)

f.write("#include <cmath>\n")
f.write("#include \"Estimation.hpp\"\n")

f.write("double polyCoeff3(double** U, double** V, double s, double r)\n{\n return %s;\n}\n\n" % cxxcode(poly.all_coeffs()[0]))
f.write("double polyCoeff2(double** U, double** V, double s, double r)\n{\n return %s;\n}\n\n" % cxxcode(poly.all_coeffs()[1]))
f.write("double polyCoeff1(double** U, double** V, double s, double r)\n{\n return %s;\n}\n\n" % cxxcode(poly.all_coeffs()[2]))
f.write("double polyCoeff0(double** U, double** V, double s, double r)\n{\n return %s;\n}\n\n" % cxxcode(poly.all_coeffs()[3]))

f.close()
