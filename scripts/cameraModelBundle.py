from sympy import *

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

t = Matrix([Symbol("Cx"), Symbol("Cy"), Symbol("Cz")])

X = Matrix([Symbol("X"), Symbol("Y"), Symbol("Z"), 1])

Extrinsics = Matrix(R).row_join(t)

P = Matrix(Extrinsics*X)

x = P[0]
y = P[1]
z = P[2]

x_ = x/z
y_ = y/z


proj = Matrix(K*Matrix([x_, y_, 1]))

u_ = proj[0]/proj[2]
v_ = proj[1]/proj[2]

f = open("src/camera/codegen/CameraModelBundle.cpp", "a")

f.truncate(0)

f.write("#include <cmath>\n")
f.write("#include \"CameraModelBundle.hpp\"\n")

f.write("double bund_u(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n" % cxxcode(u_))
f.write("double bund_v(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"% cxxcode(v_))
f.write("double bund_u_dfx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"% cxxcode(diff(u_, Symbol("fx"))))
f.write("double bund_u_dfy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("fy"))))
f.write("double bund_v_dfx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"% cxxcode(diff(v_, Symbol("fx"))))
f.write("double bund_v_dfy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("fy"))))
f.write("double bund_u_dcx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("cx"))))
f.write("double bund_u_dcy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("cy"))))
f.write("double bund_v_dcx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("cx"))))
f.write("double bund_v_dcy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("cy"))))
f.write("double bund_u_dCx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Cx"))))
f.write("double bund_u_dCy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Cy"))))
f.write("double bund_u_dCz(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Cz"))))
f.write("double bund_v_dCx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Cx"))))
f.write("double bund_v_dCy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Cy"))))
f.write("double bund_v_dCz(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Cz"))))
f.write("double bund_u_dX(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("X"))))
f.write("double bund_u_dY(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Y"))))
f.write("double bund_u_dZ(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Z"))))
f.write("double bund_v_dX(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("X"))))
f.write("double bund_v_dY(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Y"))))
f.write("double bund_v_dZ(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Z"))))
f.write("double bund_u_dk1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k1"))))
f.write("double bund_u_dk2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k2"))))
f.write("double bund_u_dk3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k3"))))
f.write("double bund_u_dk4(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k4"))))
f.write("double bund_u_dk5(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k5"))))
f.write("double bund_u_dk6(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k6"))))
f.write("double bund_u_dp1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("p1"))))
f.write("double bund_u_dp2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("p2"))))
f.write("double bund_v_dk1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("k1"))))
f.write("double bund_v_dk2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("k2"))))
f.write("double bund_v_dk3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("k3"))))
f.write("double bund_v_dp1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("p1"))))
f.write("double bund_v_dp2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("p2"))))

f.write("double bund_u_dr1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("r1"))))
f.write("double bund_u_dr2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("r2"))))
f.write("double bund_u_dr3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("r3"))))

f.write("double bund_v_dr1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("r1"))))
f.write("double bund_v_dr2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("r2"))))
f.write("double bund_v_dr3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("r3"))))

f.close()

f = open("src/camera/codegen/CameraModelBundle.hpp", "a")
f.truncate(0)
f.write("#pragma once\n")

f.write("double bund_u(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dfx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dfy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dfx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dfy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dcx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dcy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dcx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dcy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dCx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dCy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dCz(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dCx(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dCy(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dCz(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dX(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dY(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dZ(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dX(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dY(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dZ(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dk1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dk2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dk3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dk4(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dk5(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dk6(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dp1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dp2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dk1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dk2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dk3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dp1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dp2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")

f.write("double bund_u_dr1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dr2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_u_dr3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")

f.write("double bund_v_dr1(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dr2(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")
f.write("double bund_v_dr3(double fx, double fy, double cx, double cy, double Cx, double Cy, double Cz, double r1, double r2, double r3, double k1, double k2, double k3, double p1, double p2, double X, double Y, double Z);\n")

f.close()
