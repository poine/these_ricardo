#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import numpy as np, scipy.integrate
import matplotlib.pyplot as plt
import logging
import control.matlab

import single_point_load as dyn

LOG = logging.getLogger('pvtol_pole')

def sim_open_loop(P, X0, U):
    time = np.arange(0., 10, 0.01)
    X = scipy.integrate.odeint(P.dyn, X0, time, args=(U, ))
    U = U*np.ones((len(time), P.i_size))
    return time, X, U


# def sim_state_feedback(P, X0, Ue, K, H, Ycf, tf=10.):
#     def dU(X, t): return -np.dot(K, X) + np.dot(H, Ycf(t))
#     def cl_dyn(X, t): return dyn.dyn(X, t, Ue+dU(X, t), P)
#     time = np.arange(0., tf, 0.01)
#     X = scipy.integrate.odeint(cl_dyn, X0, time)
#     U =  np.array([Ue+dU(Xi, ti) for Xi, ti in zip(X, time)])
#     return time, X, U



def test_open_loop():
    P = dyn.PVTP()
    X0 = np.zeros(P.s_size)
    X0[P.s_th] = np.deg2rad(0.5 )
    X0[P.s_ph] = np.deg2rad(-1.)
    time, X, U = sim_open_loop(P, X0, P.Ue)
    dyn.plot_trajectory(time, X, U)
    anim = dyn.animate(time, X, U, P)
    plt.show()


class PvtReg:
    def __init__(self):
        pass
    
def test_state_feedback(_kind='LQR'):
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
    breakpoint()
    
    
if __name__ == "__main__":
    test_open_loop()
    #test_state_feedback()
