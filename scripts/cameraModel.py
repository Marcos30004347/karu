from sympy import *
from sympy.tensor.functions import TensorProduct

K = Matrix([[Symbol("fx"), 0,             Symbol("cx")], 
            [0,             Symbol("fy"), Symbol("cy")], 
            [0,             0,             1         ]])


r = Matrix([Symbol("r1"), Symbol("r2"), Symbol("r3")])

u = r/(sqrt(r[0]**2 + r[1]**2 + r[2]**2))

skew = Matrix([[0, -u[2], u[1]], [u[2], 0,  -u[0]], [-u[1],  u[0], 0]])

angle = sqrt(r[0]**2 + r[1]**2 + r[2]**2)

R = cos(angle)*Identity(3) + sin(angle)*skew + (1-cos(angle))*Matrix([
    [u[0]*u[0], u[0]*u[1], u[0]*u[2]],
    [u[1]*u[0], u[1]*u[1], u[1]*u[2]],
    [u[2]*u[0], u[2]*u[1], u[2]*u[2]]
])

C = Matrix([Symbol("Cx"), Symbol("Cy"), Symbol("Cz")])

X = Matrix([Symbol("X"), Symbol("Y"), Symbol("Z"), 1])

# t = Matrix(R.T*-1*C)

Extrinsics = Matrix(R.T)*Matrix(Identity(3)).row_join(-1*C)

P = Matrix(Extrinsics*X)

x = P[0]
y = P[1]
z = P[2]

x__ = x/z
y__ = y/z

# r = (x_)**2 + (y_)**2

# x__ = x_/(1 + Symbol("k1")*r + Symbol("k2")*r**2 + Symbol("k3")*r**3) + 2*Symbol("p1")*x_*y_ + Symbol("p2")*(r + 2*x_**2)
# y__ = y_/(1 + Symbol("k1")*r + Symbol("k2")*r**2 + Symbol("k3")*r**3) + 2*Symbol("p2")*x_*y_ + Symbol("p1")*(r + 2*y_**2)

# P = K*R*Matrix(Identity(3)).row_join(-1*C)

proj = Matrix(K*Matrix([x__, y__, 1]))

u_ = proj[0]/proj[2]
v_ = proj[1]/proj[2]

f = open("src/camera/codegen/cameraModel.cpp", "a")

f.truncate(0)

f.write("#include <cmath>\n")

f.write("double u(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n" % cxxcode(u_))
f.write("double v(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"% cxxcode(v_))
f.write("double u_dfx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"% cxxcode(diff(u_, Symbol("fx"))))
f.write("double u_dfy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("fy"))))
f.write("double v_dfx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"% cxxcode(diff(v_, Symbol("fx"))))
f.write("double v_dfy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("fy"))))
f.write("double u_dcx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("cx"))))
f.write("double u_dcy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("cy"))))
f.write("double v_dcx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("cx"))))
f.write("double v_dcy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("cy"))))
f.write("double u_dCx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Cx"))))
f.write("double u_dCy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Cy"))))
f.write("double u_dCz(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Cz"))))
f.write("double v_dCx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Cx"))))
f.write("double v_dCy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Cy"))))
f.write("double v_dCz(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Cz"))))
f.write("double u_dX(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("X"))))
f.write("double u_dY(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Y"))))
f.write("double u_dZ(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Z"))))
f.write("double v_dX(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("X"))))
f.write("double v_dY(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Y"))))
f.write("double v_dZ(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Z"))))
f.write("double u_dk1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k1"))))
f.write("double u_dk2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k2"))))
f.write("double u_dk3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k3"))))
f.write("double u_dk4(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k4"))))
f.write("double u_dk5(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k5"))))
f.write("double u_dk6(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k6"))))
f.write("double u_dp1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("p1"))))
f.write("double u_dp2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("p2"))))
f.write("double v_dk1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("k1"))))
f.write("double v_dk2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("k2"))))
f.write("double v_dk3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("k3"))))
f.write("double v_dp1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("p1"))))
f.write("double v_dp2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("p2"))))

f.write("double u_dr1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("r1"))))
f.write("double u_dr2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("r2"))))
f.write("double u_dr3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("r3"))))

f.write("double v_dr1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("r1"))))
f.write("double v_dr2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("r2"))))
f.write("double v_dr3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("r3"))))

f.close()

f = open("src/camera/codegen/cameraModel.hpp", "a")
f.truncate(0)

f.write("double u(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dfx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dfy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dfx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dfy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dcx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dcy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dcx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dcy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dCx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dCy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dCz(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dCx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dCy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dCz(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dX(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dY(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dZ(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dX(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dY(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dZ(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dk1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dk2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dk3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dk4(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dk5(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dk6(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dp1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dp2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dk1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dk2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dk3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dp1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dp2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")

f.write("double u_dr1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dr2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double u_dr3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")

f.write("double v_dr1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dr2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double v_dr3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")

f.close()
