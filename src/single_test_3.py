#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import sys
import numpy as np, matplotlib.pyplot as plt

import single as dyn, misc_utils as mu


import control

#
# (algebraic) Trajectory following for single vehicle
#

class CircleTraj:
    def __init__(self):
        self.c, self.r, self.om = np.array([0,0]), 1., 2.

    def get(self, t):
        a = t*self.om; ca, sa = np.cos(a), np.sin(a)
        Y0 = self.c + self.r  * np.array([ ca,  sa])
        Y1 = self.r*self.om    * np.array([-sa,  ca])
        Y2 = self.r*self.om**2 * np.array([-ca, -sa])
        Y3 = self.r*self.om**3 * np.array([ sa, -ca])
        Y4 = self.r*self.om**4 * np.array([ ca,  sa])
        return np.vstack((Y0, Y1, Y2, Y3, Y4))


class PolyTraj:
    def __init__(self):
        self.px = mu.PolynomialOne([0, 0, 0, 0, 0], [2, 0, 0, 0, 0], 10.)
    def get(self, t):
        Y = np.vstack((self.px.get(t), [0, 0, 0, 0, 0])).T
        return Y
    
def feedback_control(X, Xsp, Ue, K):
    dX = X-Xsp
    #dX[dyn.sv_th] = mu.normMpiPi(dX[dyn.sv_th])
    return Ue-K@dX

def output_to_state(Yc, P):
    Xc = np.zeros(dyn.sv_size)
    Xc[dyn.sv_slice_pos] = Yc[0]
    Xc[dyn.sv_slice_vel] = Yc[1]
    x2, z2pg = Yc[2,0], Yc[2,1]+P.g
    Xc[dyn.sv_th] = -np.arctan2(x2, z2pg)
    x3, z3 = Yc[3,0], Yc[3,1]
    Xc[dyn.sv_thd] = -(z2pg*x3-x2*z3)/(z2pg**2+x2**2)
    ut = P.m*np.sqrt(x2**2 + z2pg**2)
    x4, z4 = Yc[4,0], Yc[4,1]
    a = x4*z2pg - z4*x2
    b = x2**2 + z2pg**2
    c = 2*(z2pg*z3 + x2*x3)
    d = x3*z2pg-z3*x2
    ud = -P.J/P.d*(a/b - c*d/b**2)
    Uc = np.array([(ut+ud)/2., (ut-ud)/2.])
    return Xc, Uc

def sim(P, dt=0.01, tf=10.):
    time = np.arange(0., tf, dt)
    X, Xr, U = [np.zeros((len(time), _s)) for _s in [dyn.sv_size, dyn.sv_size, dyn.iv_size]]
    traj = CircleTraj()
    #traj = PolyTraj()
    Yc = np.array([traj.get(_t) for _t in time])
    
    Xe, Ue = dyn.trim(P)
    A, B = dyn.jacobian(Xe, Ue, P)
    Q, R = [10, 10, 1, 1, 1, 1], [1, 1]
    (K, _X, _E) = control.lqr(A, B, np.diag(Q), np.diag(R))
    K = np.array(K)
    X[0] = Xe
    for i in range(0, len(time)):
        Xr[i], Uc = output_to_state(Yc[i], P)
        U[i] =  feedback_control(X[i], Xr[i], Uc, K)
        if i < len(time)-1: X[i+1] = dyn.disc_dyn(X[i], U[i], P, dt)
    
    return time, X, U, Yc, Xr

def main(save=False):
    P = dyn.Param()  
    time, X, U, Yc, Xr = sim(P)
    dyn.plot_trajectory(time, X, U, Yc[:,0], Xr, window_title="LQR")
    if save:
        plt.savefig(mu.PLOT_DIR+f'/single_circle_tracking_chrono.png')
    anim = dyn.animate(time, X, U, P, Yc[:,0])
    if save:
        mu.save_anim(mu.PLOT_DIR+'/single_circle_tracking.apng', anim)
    plt.show()

    
if __name__ == "__main__":
    main(save='-s' in sys.argv)
