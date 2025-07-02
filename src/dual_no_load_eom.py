#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import sys, numpy as np, sympy as sp, matplotlib.pyplot as plt

import dual_no_load 

sp.init_printing(use_unicode=True)

def rmat(_a): return  sp.Matrix([[sp.cos(_a), -sp.sin(_a)],[sp.sin(_a), sp.cos(_a)]])
def n2(_a): return _a.dot(_a)

def dual_no_load_get_sym_eom_lagrange():
    t = sp.symbols('t')
    # State
    x1, z1, ph, th1, th2 = sp.symbols('x1 z1 ph th1 th2', cls=sp.Function)
    x1d, z1d, phd, th1d, th2d = [sp.diff(_e, t) for _e in [x1(t), z1(t), ph(t), th1(t), th2(t)]]
    # Input
    inputs = F1, F2, Tau1, Tau2 = sp.symbols('F1 F2 tau1 tau2')
    # Parameters
    params = m1, J1, m2, J2, l, g = sp.symbols('m1 J1 m2 J2 l g')
    
    # Master
    P1, P1d = sp.Matrix([[x1(t)],[z1(t)]]), sp.Matrix([[x1d],[z1d]])

    # Slave (Kinematic)
    L, Rph = sp.Matrix([[l],[0]]), rmat(ph(t))
    P2 = P1 + Rph@L
    P2d = sp.diff(P2, t)
    
    # Energy
    T = 0.5*m1*n2(P1d) + 0.5*J1*th1d**2 + 0.5*m2*n2(P2d) + 0.5*J2*th2d**2
    V = g*(m1*P1[1] + m2*P2[1])
    # Lagrangian
    L = T-V
    # Partial derivatives
    dL_dx1 = sp.diff(L, x1(t))
    dL_dz1 = sp.diff(L, z1(t))
    dL_dph = sp.diff(L, ph(t))
    dL_dth1 = sp.diff(L, th1(t))
    dL_dth2 = sp.diff(L, th2(t))
    dL_dx1d = sp.diff(L, sp.diff(x1(t), t))
    dL_dz1d = sp.diff(L, sp.diff(z1(t), t))
    dL_dphd = sp.diff(L, sp.diff(ph(t), t))
    dL_dth1d = sp.diff(L, sp.diff(th1(t), t))
    dL_dth2d = sp.diff(L, sp.diff(th2(t), t))

    
    # Forces and Moments
    Fx = -F1*sp.sin(th1(t)) - F2*sp.sin(th2(t))
    Fz =  F1*sp.cos(th1(t)) + F2*sp.cos(th2(t))
    Mth, Mth2 = Tau1, Tau2
    Mph = F2*l*sp.cos(ph(t)-th2(t))
    # Lagrange
    La1 = sp.diff(dL_dx1d, t) - dL_dx1 - Fx
    x1dd1 = sp.solve(La1, sp.Derivative(x1(t), (t, 2)))[0]   # solve La1 for x1dd
    La2 = sp.diff(dL_dz1d, t) - dL_dz1 - Fz 
    z1dd1 = sp.solve(La2, sp.Derivative(z1(t), (t, 2)))[0]   # solve La2 for z1dd
    La3 = sp.diff(dL_dth1d, t) - dL_dth1 - Mth
    th1dd1 = sp.solve(La3, sp.Derivative(th1(t), (t, 2)))[0] # solve La3 for th1dd
    La4 = sp.diff(dL_dth2d, t) - dL_dth2 - Mth2
    th2dd1 = sp.solve(La4, sp.Derivative(th2(t), (t, 2)))[0] # solve La4 for th2dd
    La5 = sp.diff(dL_dphd, t) - dL_dph - Mph
    phdd1 = sp.solve(La5, sp.Derivative(ph(t), (t, 2)))[0]   # solve La5 for phdd 
    phdd2 = phdd1.subs(sp.Derivative(x1(t), (t, 2)), x1dd1).subs(sp.Derivative(z1(t), (t, 2)), z1dd1).simplify()
    phdd3 = sp.solve(phdd2, sp.Derivative(ph(t), (t, 2)))[0] # final phdd
    # subs phdd in xdd and zdd expressions
    x1dd2 = x1dd1.subs(sp.Derivative(ph(t), (t, 2)), phdd3)
    z1dd2 = z1dd1.subs(sp.Derivative(ph(t), (t, 2)), phdd3)

    # State vector x, z, th, ph, th2, xd, zd, thd, phd, th2d
    states = sp.Matrix([x1(t), z1(t), ph(t), th1(t), th2(t),
                        x1d, z1d, phd, th1d, th2d])
    eom = sp.Matrix([x1d, z1d, phd, th1d, th2d,
                     x1dd2, z1dd2, phdd3, th1dd1, th2dd1])


    return states, inputs, params, eom
   

def compare(pvt, X, U, states, inputs, params, eom):
    Xd = np.expand_dims(pvt.dyn2(X, 0., U), 1)
    # substitube parameters numerical values
    for _p, _v in zip(params, [pvt.P.m1, pvt.P.J1, pvt.P.m2, pvt.P.J2, pvt.P.l, pvt.P.g]):
        eom = eom.subs(_p, _v)
    # substitute input numerical values
    F1, F2 = U[0]+U[1], U[2]+U[3]
    Tau1, Tau2 = (-U[0]+U[1])*pvt.P.d1, (-U[2]+U[3])*pvt.P.d2 
    for _i, _v in zip(inputs, [F1, F2, Tau1, Tau2]):
        eom = eom.subs(_i, _v)
    # substitue state
    for _i in range(5, pvt.s_size):
        eom = eom.subs(states[_i], X[_i])
    for _i in range(5):
        eom = eom.subs(states[_i], X[_i])
    print(np.allclose(Xd, np.array(eom).astype(np.float64)))
    return Xd, eom
    
def main():
    states, inputs, params, eom = dual_no_load_get_sym_eom_lagrange()
    #breakpoint()
    #print(eom)

    #eom_sym = eom
    pvt = dual_no_load.PVT()

    X = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    U = np.array(pvt.Ue)
    
    X = [0, 0, 0, 0, 0, 0.1, 0.2, 0.1, 0, 0]
    U = np.array(pvt.Ue)
    #U[pvt.i_fl1] *= 1.5
    #U[pvt.i_fr1] *= 1.5
    U[pvt.i_fl1] *= 1.1
    U[pvt.i_fl2] *= 1.1
    
    Xd, eom = compare(pvt, X, U, states, inputs, params, eom)

    
    #print(Xd)
    #print(repr(eom))
    #breakpoint()

    
if __name__ == "__main__":
    main()

