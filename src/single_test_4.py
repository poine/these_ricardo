#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import sys
import numpy as np, matplotlib.pyplot as plt
import control

import single as dyn, misc_utils as mu
import single_test_3 as test_3

#
# Reference Model Trajectory following for single vehicle
#

def stepxz(t):
    Yc = np.zeros(2)
    #Yc = mu.step(t, 5.), -mu.step(t, 2.5)
    Yc = mu.step(t, 5.), 0
    return Yc

def feedback_control(X, Xsp, Ur, K):
    dX = X-Xsp
    #dX[PVT.PVT.s_th] = normMpiPi(dX[PVT.PVT.s_th])
    return -K@dX


class TwoDRef:
    def __init__(self):
        poles = [complex(-3, 3), complex(-3, -3), complex(-4, 4), complex(-4, -4), -5]
        coefs = np.poly(poles)
        _sats = [5., 10., 100., 500, 15000.]  # vel, accel, jerk, snap, crackle 
        self.refx = mu.LinRef(-np.flip(coefs)[:-1], sats=_sats)
        self.refz = mu.LinRef(-np.flip(coefs)[:-1], sats=_sats)

    def get(self, sp, tim1, ti):
        Yr = np.zeros((6,2))
        Yr[:,0] = self.refx.run(ti-tim1, sp[0])
        Yr[:,1] = self.refz.run(ti-tim1, sp[1])
        return Yr
        
def sim(P, dt=0.01, tf=10.):
    time = np.arange(0., tf, dt)
    X, Xr, Ur, U = [np.zeros((len(time), _s)) for _s in [dyn.sv_size, dyn.sv_size, dyn.iv_size, dyn.iv_size]]
    Yr = np.zeros((len(time), 6, 2))
    r = TwoDRef()
    Xe, Ue = dyn.trim(P)
    A, B = dyn.jacobian(Xe, Ue, P)
    Q, R = [5, 10, 1, 1, 1, 1], [1, 1]
    (K, _X, _E) = control.lqr(A, B, np.diag(Q), np.diag(R))
    K = np.array(K)
    X[0] = Xr[0]
    for i in range(0, len(time)):
        if i > 1:
            Yr[i] = r.get(stepxz(time[i]), time[i-1], time[i])
        Xr[i], Ur[i] = test_3.output_to_state(Yr[i], P)
        #breakpoint()
        U[i] = np.flip(Ur[i]) + feedback_control(X[i], Xr[i], Ur[i], K)
        if i < len(time)-1:
            X[i+1] = dyn.disc_dyn(X[i], U[i], P, dt)
        
    dyn.plot_trajectory(time, X, Ur, Yr[:,0], Xr, window_title="LQR")
    anim = dyn.animate(time, X, U, Yr[:,0], P)
    plt.show()
    
def main(save=False):
    P = dyn.Param()  
    sim(P)
    
if __name__ == "__main__":
    main(save='-s' in sys.argv)

    
