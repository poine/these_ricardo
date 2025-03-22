#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#
#
# Dual tethered, master-slave feedback regulation
#
#

import logging, sys
import numpy as np, scipy.integrate, matplotlib.pyplot as plt
import control.matlab

import dual_no_load as PVT, misc_utils as mu

LOG = logging.getLogger(__name__)


def ms_feedback(X, Xsp, Ue, K1, K2):
    U = np.array(Ue)
    dX_master = PVT.PVT.master_state(X-Xsp)
    dU1 = - K1@dX_master # master vehicle
    dX_slave = PVT.PVT.slave_state(X-Xsp)
    dU2 = K2@dX_slave
    #breakpoint()
    return Ue-dU1-dU2

def feedback_control(X, Xsp, Ue, K):
    return Ue - K@(X-Xsp)

def stepx(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_x] = mu.step(t, 0.5)
    return X
def cst(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_x] = -0.25
    return X

def sim_feedback(sp, dt=0.01):
    P = PVT.PVT()
    X0, Ue = np.zeros(P.s_size), P.Ue
    if 0:
        # feedbak regulation for master drone
        import single
        _p = single.Param()
        _A, _B = single.jacobian(*single.trim(_p), _p)
        _Q, _R = [3, 10, 0.5, 0.25, 0.25, 0.005], [1, 1]
        (_K1, X, E) = control.lqr(_A, _B, np.diag(_Q), np.diag(_R))
        print(f"master's poles {np.linalg.eig(_A-_B@_K1)[0]}")
        K1 = np.vstack((_K1, np.zeros((2,6))))
        _tmp = np.array(K1[0,:]) # single is defined backward FIXME
        K1[0,:] = K1[1,:]
        K1[1,:] = _tmp 
        #breakpoint()
        K2 = np.zeros((4,4))
        #K2[2,:] = [2.,  -0.0, 0.25,  -0.0]
        #K2[3,:] = [2.,   0.0, 0.25,   0.0]
        K2[2,:] = [.01,   0.0, 0.08,   0.0]
        K2[3,:] = [.01,   0.0, 0.08,   0.0]
    else:
        A, B = P.jac()
        Q = np.diag([3., 7.,0.5, 10., 0.5,
                     0.25, 0.25, 0.005, 0.25, 0.005])
        R = np.diag([0.5, 0.5, 0.5, 0.5])
        (K, X, E) = control.matlab.lqr(A, B, Q, R)
        poles, vect_p = np.linalg.eig(A-np.dot(B, K))
        with np.printoptions(precision=2, linewidth=200):
            LOG.info('gain:\n{}'.format(K))
            LOG.info('poles:\n{}'.format(poles))
    
    X0 =np.array([0.1, 0.2, np.deg2rad(1),     np.deg2rad(1), np.deg2rad(2),
                  0.01, 0.02, np.deg2rad(0.1), np.deg2rad(3), np.deg2rad(4)])
    X0 =np.zeros(10)
    X0[PVT.PVT.s_x] = -0.25
    X0[PVT.PVT.s_ph] = np.deg2rad(10)
    time = np.arange(0., 4., dt)
    X = np.zeros((len(time), P.s_size))
    Xsp = np.zeros((len(time), P.s_size))
    U = np.zeros((len(time), P.i_size))
    X[0] = X0
    for i in range(0, len(time)):
        Xsp[i] = sp(time[i])
        #U[i] =  ms_feedback(X[i], Xsp[i], Ue, K1, K2)
        U[i] =  feedback_control(X[i], Xsp[i], Ue, K)
        if i < len(time)-1: X[i+1] = P.disc_dyn(X[i], U[i], dt)
    PVT.plot_trajectory(time, X, U)
    anim = PVT.animate(time, X, U, P)
    plt.show()
    
def main(exp):
    logging.basicConfig(level=logging.INFO)
    sim_feedback(cst)
    #sim_feedback(stepx)
    
if __name__ == "__main__":
    exp = 0 if len(sys.argv)<2 else int(sys.argv[1])
    main(exp)
