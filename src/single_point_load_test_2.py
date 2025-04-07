#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import numpy as np, scipy.integrate
import matplotlib.pyplot as plt
import logging, sys
import control.matlab

import single_point_load as dyn, misc_utils as mu

LOG = logging.getLogger('pvtol_pole')

def stepx(t):
    Xsp = np.zeros(dyn.PVTP.s_size)
    Xsp[dyn.PVTP.s_x] = mu.step(t, 1.)
    return Xsp

def state_feedback(X, Xe, Ue, K): return Ue - K@(X-Xe)

def sim_state_feedback(X0, time, sp, Ue, K, vehicule):
    dt, nt = time[1]-time[0], len(time)
    X, Xsp, U = [np.zeros((nt, _s)) for _s in [vehicule.s_size, vehicule.s_size, vehicule.i_size]] 
    X[0] = X0
    for i in range(len(time)):
        Xsp[i] = sp(time[i])
        U[i] = state_feedback(X[i], Xsp[i], Ue, K)
        if i<len(time)-1:
            X[i+1] = vehicule.disc_dyn(X[i], U[i], dt)
    return time, X, Xsp, U
    
def test_state_feedback(save=False, sf='single_point_load__state_feedback_1', _kind='LQR'):
    P = dyn.PVTP()
    Xe = np.zeros(P.s_size)
    Ue = P.P.mt*P.P.g/2*np.ones(2)
    A, B = P.jac()
    if _kind == 'place':
        poles=[-5, -5, -5, -5, -5, -5, -5, -5]
        self.K = control.matlab.place(A, B, poles)
    else:
        Q = np.diag([1., 10., 1.5, 0.25, 0.01, 0.05, 0.005, 0.001])
        R = np.diag([0.5, 0.5])
        (K, X, E) = control.matlab.lqr(A, B, Q, R)
    poles, vect_p = np.linalg.eig(A-np.dot(B, K))
    LOG.info('gain:\n{}'.format(K))
    LOG.info('poles:\n{}'.format(poles))
    time, X0 = np.arange(0, 10, 0.01), np.array(Xe)
    X0[P.s_x] += 1.
    time, X, Xsp, U = sim_state_feedback(X0, time, stepx, Ue, K, P)
    dyn.plot_trajectory(time, X, U)
    if save:
        plt.savefig(mu.PLOT_DIR+f'/{sf}_chrono.png')
        anim = dyn.animate(time, X, U, P)
    if save:
        mu.save_anim(mu.PLOT_DIR+f'/{sf}.apng', anim)
    plt.show()
    #breakpoint()
    
    
if __name__ == "__main__":
    test_state_feedback(save='-s' in sys.argv)
