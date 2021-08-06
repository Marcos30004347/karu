from sympy import *

U = Matrix([
    [Symbol("U1"), 0, 0],
    [0, Symbol("U2"), 0],
    [0, 0, Symbol("U3")],
])

W = Matrix([
    [Symbol("W11"), Symbol("W12"), Symbol("W13"), Symbol("W14")],
    [Symbol("W21"), Symbol("W22"), Symbol("W23"), Symbol("W24")],
    [Symbol("W31"), Symbol("W32"), Symbol("W33"), Symbol("W34")],
])

V_inv = Matrix([
    [Symbol("V1^-1"), 0, 0, 0],
    [0, Symbol("V2^-1"), 0, 0],
    [0, 0, Symbol("V3^-1"), 0],
    [0, 0, 0, Symbol("V4^-1")],
])

W_T = Matrix([
    [transpose(Symbol("W11")), transpose(Symbol("W21")), transpose(Symbol("W31"))],
    [transpose(Symbol("W12")), transpose(Symbol("W22")), transpose(Symbol("W32"))],
    [transpose(Symbol("W13")), transpose(Symbol("W23")), transpose(Symbol("W33"))],
    [transpose(Symbol("W14")), transpose(Symbol("W24")), transpose(Symbol("W34"))],
])

init_printing(wrap_line=False)

pprint(U - W*V_inv*W_T)

rc = Matrix([Symbol("rc1"), Symbol("rc2"), Symbol("rc3")])
rp = Matrix([Symbol("rp1"), Symbol("rp2"), Symbol("rp3"), Symbol("rp4")])

pprint(rc - W*V_inv*rp)
