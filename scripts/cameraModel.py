from sympy import *
import numpy as np
import math

# This function return an rotation matrix from and axes vector and angle
def get_rotation_matrix(axis, angle):
    rad = math.radians(angle)
    vec = rad * axis/np.linalg.norm(axis)

    r = Matrix([Symbol("r1"), Symbol("r2"), Symbol("r3")])
    a = Symbol("a")

    skew_cross = Matrix([[0, -r[2]/a, r[1]/a], [r[2]/a, 0,  -r[0]/a], [-r[1]/a,  r[0]/a, 0]])

    R = cos(a)*Identity(3) + sin(a)*skew_cross + (1-cos(a))*((r/a)*(r.T/a))

    rotation = Matrix(R).expand()

    rotation = rotation.subs(a, rad)
    rotation = rotation.subs(r[0], vec[0])
    rotation = rotation.subs(r[1], vec[1])
    rotation = rotation.subs(r[2], vec[2])

    return rotation

def camera_model():
  # angle of the rotation on the axis r
  # a = Symbol("a")

  K = Matrix([[Symbol("fx"), 0,             Symbol("cx")], 
              [0,             Symbol("fy"), Symbol("cy")], 
              [0,             0,             1            ]])

  C = Matrix([Symbol("Cx"), Symbol("Cy"), Symbol("Cz")])

  # Matrix containing the axis and angle of the rotation, axis = r/|r| and anngle = |r|
  # the angle also needs to be equal to the parameter a, so a = |r|
  r = Matrix([Symbol("r1"), Symbol("r2"), Symbol("r3")])

  # Skew-Symetrix-Matrix for the cross produc of r
  rc = Matrix([[0, -r[2], r[1]], [r[2], 0,  -r[0]], [-r[1],  r[0], 0]])
  rc = rc/(sqrt(r[0]**2 + r[1]**2 + r[2]**2))
  # Rotation matrix
  R = cos(sqrt(r[0]**2 + r[1]**2 + r[2]**2))*Identity(3) + sin(sqrt(r[0]**2 + r[1]**2 + r[2]**2))*rc + (1-cos(sqrt(r[0]**2 + r[1]**2 + r[2]**2)))*((r/sqrt(r[0]**2 + r[1]**2 + r[2]**2))*(r.T/sqrt(r[0]**2 + r[1]**2 + r[2]**2)))
  
  # T = Matrix(Identity(3)).row_join(Matrix(-1*C))

  # P = Matrix([Symbol("X"), Symbol("Y"), Symbol("Z"), 1])
  
  # Proj = K*R*T*P

  # projection = Matrix(Proj).expand()

  # return projection
  # projection = Matrix(Proj).expand()

  P = Matrix([Symbol("X"), Symbol("Y"), Symbol("Z")])
  
  v = Matrix((R*P + C))

  x = v[0]
  y = v[1]
  z = v[2]

  x_ = x/z
  y_ = y/z

  r = (x_)**2 + (y_)**2

  x__ = x_ * (1 + Symbol("k1")*r + Symbol("k2")*r**2 + Symbol("k3")*r**3) + 2*Symbol("p1")*x_*y_ + Symbol("p2")*(r + 2*x_**2)
  y__ = y_ * (1 + Symbol("k1")*r + Symbol("k2")*r**2 + Symbol("k3")*r**3) + 2*Symbol("p2")*x_*y_ + Symbol("p1")*(r + 2*y_**2)

  View = K*Matrix([x__, y__, 1])

  return (View[0], View[1])

u, v = camera_model()

u_ = u
v_ = v

print("#pragma once")

print("#include <cmath>")

print("// camera focus parameters derivatives")
print("float u(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return", cxxcode(u_), end=";\n}\n\n")
print("float u_dfx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return", cxxcode(diff(u_, Symbol("fx"))), end=";\n}\n\n")
print("float u_dfy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("fy"))), end=";\n}\n\n")
print("float v(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return", cxxcode(v_), end=";\n}\n\n")
print("float v_dfx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return", cxxcode(diff(v_, Symbol("fx"))), end=";\n}\n\n")
print("float v_dfy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("fy"))), end=";\n}\n\n")

print("// camera principal point derivatives")
print("float u_dcx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("cx"))), end=";\n}\n\n")
print("float u_dcy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("cy"))), end=";\n}\n\n")
print("float v_dcx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("cx"))), end=";\n}\n\n")
print("float v_dcy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("cy"))), end=";\n}\n\n")

print("// camera position parameters derivatives")
print("float u_dCx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("Cx"))), end=";\n}\n\n")
print("float u_dCy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("Cy"))), end=";\n}\n\n")
print("float u_dCz(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("Cz"))), end=";\n}\n\n")
print("float v_dCx(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("Cx"))), end=";\n}\n\n")
print("float v_dCy(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("Cy"))), end=";\n}\n\n")
print("float v_dCz(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("Cz"))), end=";\n}\n\n")

print("// camera rotation parameters derivatives")
print("float u_da(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("a"))), end=";\n}\n\n")
print("float u_dr1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("r1"))), end=";\n}\n\n")
print("float u_dr2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("r2"))), end=";\n}\n\n")
print("float u_dr3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("r3"))), end=";\n}\n\n")
print("float v_da(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("a"))), end=";\n}\n\n")
print("float v_dr1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("r1"))), end=";\n}\n\n")
print("float v_dr2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("r2"))), end=";\n}\n\n")
print("float v_dr3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("r3"))), end=";\n}\n\n")

print("// Point parameters derivatives")
print("float u_dX(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("X"))), end=";\n}\n\n")
print("float u_dY(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("Y"))), end=";\n}\n\n")
print("float u_dZ(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("Z"))), end=";\n}\n\n")
print("float v_dX(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("X"))), end=";\n}\n\n")
print("float v_dY(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("Y"))), end=";\n}\n\n")
print("float v_dZ(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("Z"))), end=";\n}\n\n")

print("// distortion parameters derivatives")
print("// for pinhole cameras, k1=0,k2=0,k3=0,k4=0,k5=0,k6=0,p1=0,p2=0")
print("float u_dk1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("k1"))), end=";\n}\n\n")
print("float u_dk2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("k2"))), end=";\n}\n\n")
print("float u_dk3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("k3"))), end=";\n}\n\n")
print("float u_dk4(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("k4"))), end=";\n}\n\n")
print("float u_dk5(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("k5"))), end=";\n}\n\n")
print("float u_dk6(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("k6"))), end=";\n}\n\n")
print("float u_dp1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("p1"))), end=";\n}\n\n")
print("float u_dp2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(u_, Symbol("p2"))), end=";\n}\n\n")
print("float v_dk1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("k1"))), end=";\n}\n\n")
print("float v_dk2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("k2"))), end=";\n}\n\n")
print("float v_dk3(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("k3"))), end=";\n}\n\n")
print("float v_dp1(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("p1"))), end=";\n}\n\n")
print("float v_dp2(float fx, float fy, float cx, float cy, float Cx, float Cy, float Cz, float r1, float r2, float r3, float k1, float k2, float k3, float p1, float p2, float X, float Y, float Z)\n{\n return",cxxcode(diff(v_, Symbol("p2"))), end=";\n}\n\n")



# print(camera_model()[2])
# print()


# u3 = (Symbol("X") - Symbol("Cx")) * ((Symbol("r1")*Symbol("r3")*(1-cos(Symbol("a"))))/(Symbol("a")*Symbol("a")) - ((Symbol("r2")*sin(Symbol("a")))/Symbol("a"))) + (Symbol("y") - Symbol("Cy"))*(((Symbol("r2")*Symbol("r3")*(1-cos(Symbol("a"))))/(Symbol("a")*Symbol("a"))) + ((Symbol("r1")*sin(Symbol("a")))/Symbol("a"))) + (Symbol("z") - Symbol("Cz"))*( (((Symbol("r3")*Symbol("r3"))*(1 - cos(Symbol("a"))))/(Symbol("a")*Symbol("a"))) + cos(Symbol("a")))                                                                          
# print(u3.expand())
