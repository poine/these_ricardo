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
    Yc = mu.step(t, 5.), -mu.step(t, 1., p=5.)
    #Yc = mu.step(t, 5.), 0
    return Yc

def feedback_control(X, Xsp, Ur, K):
    dX = X-Xsp
    #dX[PVT.PVT.s_th] = normMpiPi(dX[PVT.PVT.s_th])
    return -K@dX


class TwoDRef:
    def __init__(self):
        poles = [complex(-3, 2), complex(-3, -2), complex(-4, 3), complex(-4, -3), -6]
        coefs = np.poly(poles)
        _sats = [5., 10., 100., 500, 15000.]  # vel, accel, jerk, snap, crackle 
        self.refx = mu.LinRef(-np.flip(coefs)[:-1], sats=_sats)
        self.refz = mu.LinRef(-np.flip(coefs)[:-1], sats=_sats)

    def get(self, sp, tim1, ti):
        Yr = np.zeros((6,2))
        Yr[:,0] = self.refx.run(ti-tim1, sp[0])
        Yr[:,1] = self.refz.run(ti-tim1, sp[1])
        return Yr
        
def sim(P, dt=0.01, tf=12., save=None):
    time = np.arange(0., tf, dt)
    X, Xr, Ur, U = [np.zeros((len(time), _s)) for _s in [dyn.sv_size, dyn.sv_size, dyn.iv_size, dyn.iv_size]]
    Ysp, Yr = np.zeros((len(time), 2)), np.zeros((len(time), 6, 2))
    r = TwoDRef()
    Xe, Ue = dyn.trim(P)
    A, B = dyn.jacobian(Xe, Ue, P)
    Q, R = [5, 10, 1, 1, 1, 1], [1, 1]
    (K, _X, _E) = control.lqr(A, B, np.diag(Q), np.diag(R))
    K = np.array(K)
    X[0] = Xr[0]
    for i in range(0, len(time)):
        if i > 1:
            Ysp[i] = stepxz(time[i])
            Yr[i] = r.get(Ysp[i], time[i-1], time[i])
        Xr[i], Ur[i] = test_3.output_to_state(Yr[i], P)
        #breakpoint()
        U[i] = np.flip(Ur[i]) + feedback_control(X[i], Xr[i], Ur[i], K)
        if i < len(time)-1:
            X[i+1] = dyn.disc_dyn(X[i], U[i], P, dt)
        
    dyn.plot_trajectory(time, X, Ur, Yr[:,0], Ysp, Xr, window_title="LQR")
    if save:
        plt.savefig(mu.PLOT_DIR+f'/single_ref_model_{save}_chrono.png')
    anim = dyn.animate(time, X, U, P, Ysp, Yr[:,0])
    if save:
        mu.save_anim(mu.PLOT_DIR+f'/single_ref_model_{save}.apng', anim)
    plt.show()
    
def main(save=False):
    P = dyn.Param()  
    sim(P, save='step_xz' if save else None)
    
if __name__ == "__main__":
    main(save='-s' in sys.argv)

    
