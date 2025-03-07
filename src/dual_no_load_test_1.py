#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#
#
# Dual tethered, open loop
#
#


import numpy as np, scipy.integrate, matplotlib.pyplot as plt


import dual_no_load as PVT, misc_utils as mu

def stepUt(P, t):
    U = np.array(P.Ue)
    ud1 = 0.005*np.sin(t+np.pi/2)
    ut1 = 0.01
    U[P.i_fl1] += ud1+ut1; U[P.i_fr1] += -ud1+ut1
    ud2, ut2 = -ud1, ut1
    U[P.i_fl2] += ud2+ut1; U[P.i_fr2] += -ud2+ut2
    return U

def test_open_loop(dt=0.01):
    P = PVT.PVT()
    X0, Ue = np.zeros(P.s_size), P.Ue
    X0[P.s_th] = np.deg2rad(11.5)
    X0[P.s_th2] = np.deg2rad(-11.5)
    time = np.arange(0., 12.7, dt)
    X = np.zeros((len(time), P.s_size))
    U = Ue*np.ones((len(time), P.i_size))
    X[0] = X0
    for i in range(0, len(time)-1):
        U[i] =  stepUt(P, time[i])
        X[i+1] = P.disc_dyn(X[i], U[i], dt)
    U[-1] =  stepUt(P, time[-1])
    PVT.plot_trajectory(time, X, U)
    anim = PVT.animate(time, X, U, P)
    plt.show()
    
def main():
    test_open_loop()


if __name__ == "__main__":
    main()
