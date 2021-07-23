from sympy import *

K = Matrix([[Symbol("fx"), 0,                Symbol("cx")], 
            [0,             -1*Symbol("fy"), Symbol("cy")], 
            [0,             0,             1            ]])

C = Matrix([Symbol("Cx"), Symbol("Cy"), Symbol("Cz")])

r = Matrix([Symbol("r1"), Symbol("r2"), Symbol("r3")])
rc = Matrix([[0, -r[2], r[1]], [r[2], 0,  -r[0]], [-r[1],  r[0], 0]])
rc = rc/(sqrt(r[0]**2 + r[1]**2 + r[2]**2))
R = cos(sqrt(r[0]**2 + r[1]**2 + r[2]**2))*Identity(3) + sin(sqrt(r[0]**2 + r[1]**2 + r[2]**2))*rc + (1-cos(sqrt(r[0]**2 + r[1]**2 + r[2]**2)))*((r/sqrt(r[0]**2 + r[1]**2 + r[2]**2))*(r.T/sqrt(r[0]**2 + r[1]**2 + r[2]**2)))

V = Matrix(R).row_join(Matrix(C))

P = Matrix([Symbol("X"), Symbol("Y"), Symbol("Z"), 1])

P = Matrix(V*P)

x = P[0]
y = P[1]
z = P[2]

x_ = x/z
y_ = y/z

r = (x_)**2 + (y_)**2

x__ = x_/(1 + Symbol("k1")*r + Symbol("k2")*r**2 + Symbol("k3")*r**3) + 2*Symbol("p1")*x_*y_ + Symbol("p2")*(r + 2*x_**2)
y__ = y_/(1 + Symbol("k1")*r + Symbol("k2")*r**2 + Symbol("k3")*r**3) + 2*Symbol("p2")*x_*y_ + Symbol("p1")*(r + 2*y_**2)

View = K*Matrix([x__, y__, 1])

u_ = View[0]
v_ = View[1]

f = open("src/camera/codegen/cameraModel.cpp", "a")
f.truncate(0)

f.write("#include <cmath>\n")

f.write("float u(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n" % cxxcode(u_))
f.write("float u_dfx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"% cxxcode(diff(u_, Symbol("fx"))))
f.write("float u_dfy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("fy"))))
f.write("float v(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"% cxxcode(v_))
f.write("float v_dfx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"% cxxcode(diff(v_, Symbol("fx"))))
f.write("float v_dfy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("fy"))))
f.write("float u_dcx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("cx"))))
f.write("float u_dcy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("cy"))))
f.write("float v_dcx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("cx"))))
f.write("float v_dcy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("cy"))))
f.write("float u_dCx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Cx"))))
f.write("float u_dCy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Cy"))))
f.write("float u_dCz(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Cz"))))
f.write("float v_dCx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Cx"))))
f.write("float v_dCy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Cy"))))
f.write("float v_dCz(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Cz"))))
f.write("float u_da(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("a"))))
f.write("float u_dr1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("r1"))))
f.write("float u_dr2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("r2"))))
f.write("float u_dr3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("r3"))))
f.write("float v_da(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("a"))))
f.write("float v_dr1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("r1"))))
f.write("float v_dr2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("r2"))))
f.write("float v_dr3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("r3"))))
f.write("float u_dX(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("X"))))
f.write("float u_dY(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Y"))))
f.write("float u_dZ(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("Z"))))
f.write("float v_dX(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("X"))))
f.write("float v_dY(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Y"))))
f.write("float v_dZ(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("Z"))))
f.write("float u_dk1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k1"))))
f.write("float u_dk2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k2"))))
f.write("float u_dk3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k3"))))
f.write("float u_dk4(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k4"))))
f.write("float u_dk5(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k5"))))
f.write("float u_dk6(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("k6"))))
f.write("float u_dp1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("p1"))))
f.write("float u_dp2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(u_, Symbol("p2"))))
f.write("float v_dk1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("k1"))))
f.write("float v_dk2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("k2"))))
f.write("float v_dk3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("k3"))))
f.write("float v_dp1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("p1"))))
f.write("float v_dp2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return %s;\n}\n\n"%cxxcode(diff(v_, Symbol("p2"))))

f.close()

f = open("src/camera/codegen/cameraModel.hpp", "a")
f.truncate(0)

f.write("float u(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dfx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dfy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dfx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dfy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dcx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dcy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dcx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dcy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dCx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dCy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dCz(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dCx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dCy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dCz(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_da(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dr1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dr2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dr3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_da(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dr1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dr2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dr3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dX(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dY(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dZ(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dX(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dY(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dZ(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dk1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dk2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dk3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dk4(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dk5(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dk6(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dp1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float u_dp2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dk1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dk2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dk3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dp1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")
f.write("float v_dp2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z);\n")

f.close()
