from typing import final
from sympy import *
import numpy as np
import math

r1 = Symbol("r1")
r2 = Symbol("r2")
r3 = Symbol("r3")

a = Symbol("a")

r = Matrix([r1, r2, r3])

a_skew_cross = Matrix([[0, -r[2]/a, -r[1]/a],
                       [r[2]/a, 0,  -r[0]/a],
                       [r[1]/a,  r[0]/a, 0]])

R = cos(a)*Identity(3) + sin(a)*a_skew_cross + (1-cos(a))*((r/a)*(r.T/a))
rotation = Matrix(R).expand()

rotation = rotation.subs(a, np.sqrt(2))
rotation = rotation.subs(r1, 1)
rotation = rotation.subs(r2, 0)
rotation = rotation.subs(r3, 1)

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

#print(get_rotation_matrix(np.array([0.14, 0.5, 0.75]), 45))

def camera_model():
    # focal length
    f = Symbol("f")

    # angle of the rotation on the axis r
    a = Symbol("a")

    K = Matrix([[Symbol("k11"), 0,             Symbol("k13")], 
                [0,             Symbol("k11"), Symbol("k23")], 
                [0,             0,             1            ]])

    C = Matrix([Symbol("Cx"), Symbol("Cy"), Symbol("Cz")])

    # Matrix containing the axis and angle of the rotation, axis = r/|r| and anngle = |r|
    # the angle also needs to be equal to the parameter a, so a = |r|
    r = Matrix([Symbol("r1"), Symbol("r2"), Symbol("r3")])

    # Skew-Symetrix-Matrix for the cross produc of r
    rc = Matrix([[0, -r[2]/a, r[1]/a], [r[2]/a, 0,  -r[0]/a], [-r[1]/a,  r[0]/a, 0]])

    # Rotation matrix
    R = cos(a)*Identity(3) + sin(a)*rc + (1-cos(a))*((r/a)*(r.T/a))
    T = Matrix(Identity(3))
    T = T.row_join(Matrix(-1*C))

    P = Matrix([Symbol("X"), Symbol("Y"), Symbol("Z"), 1])
    Proj = K*R*T*P
    projection = Matrix(Proj).expand()

    return projection

print(camera_model()[2])
print()
#print(simplify(camera_model()[2]))


u3 = (Symbol("X") - Symbol("Cx")) * ((Symbol("r1")*Symbol("r3")*(1-cos(Symbol("a"))))/(Symbol("a")*Symbol("a")) - ((Symbol("r2")*sin(Symbol("a")))/Symbol("a"))) + (Symbol("y") - Symbol("Cy"))*(((Symbol("r2")*Symbol("r3")*(1-cos(Symbol("a"))))/(Symbol("a")*Symbol("a"))) + ((Symbol("r1")*sin(Symbol("a")))/Symbol("a"))) + (Symbol("z") - Symbol("Cz"))*( (((Symbol("r3")*Symbol("r3"))*(1 - cos(Symbol("a"))))/(Symbol("a")*Symbol("a"))) + cos(Symbol("a")))                                                                          
print(u3.expand())
#(a**2*(-Cz + Z)*cos(a) + a*(Cx*r2 - Cy*r1 - X*r2 + Y*r1)*sin(a) + r3*(Cx*r1*cos(a) - Cx*r1 + Cy*r2*cos(a) - Cy*r2 + Cz*r3*cos(a) - Cz*r3 - X*r1*cos(a) + X*r1 - Y*r2*cos(a) + Y*r2 - Z*r3*cos(a) + Z*r3))/a**2
