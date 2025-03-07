#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#
#
# Dual tethered, feedback regulation
#
#

import logging, sys
import numpy as np, scipy.integrate, matplotlib.pyplot as plt
import control.matlab

import dual_no_load as PVT, misc_utils as mu

LOG = logging.getLogger(__name__)

def feedback_control(X, Xsp, Ue, K):
    return Ue - K@(X-Xsp)

def stepx(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_x] = mu.step(t, 0.5)
    return X

def stepz(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_z] = mu.step(t, 0.4)
    return X

def stepphi(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_ph] = mu.step(t, np.deg2rad(15))
    return X


def sim_feedback(sp, dt=0.01):
    P = PVT.PVT()
    X0, Ue = np.zeros(P.s_size), P.Ue
    A, B = P.jac()
    #    s_x, s_z, s_ph, s_th, s_th2
    Q = np.diag([3., 7., 10., 0.5, 0.5,
                 0.25, 0.25, 0.25, 0.005, 0.005])
    R = np.diag([0.5, 0.5, 0.5, 0.5])
    (K, X, E) = control.matlab.lqr(A, B, Q, R)
    poles, vect_p = np.linalg.eig(A-np.dot(B, K))
    LOG.info('gain:\n{}'.format(K))
    LOG.info('poles:\n{}'.format(poles))

    X0[P.s_th] = np.deg2rad(11.5)
    X0[P.s_th2] = np.deg2rad(-11.5)
    time = np.arange(0., 12.7, dt)
    X = np.zeros((len(time), P.s_size))
    Xsp = np.zeros((len(time), P.s_size))
    U = np.zeros((len(time), P.i_size))
    X[0] = X0
    for i in range(0, len(time)):
        Xsp[i] = sp(time[i])
        U[i] =  feedback_control(X[i], Xsp[i], Ue, K)
        if i < len(time)-1: X[i+1] = P.disc_dyn(X[i], U[i], dt)
    PVT.plot_trajectory(time, X, U)
    anim = PVT.animate(time, X, U, P)
    plt.show()


def test_1(): sim_feedback(stepx)
def test_2(): sim_feedback(stepz)
def test_3(): sim_feedback(stepphi)

def main(exp):
    logging.basicConfig(level=logging.INFO)
    exps = [test_1, test_2, test_3]
    exps[exp]()

if __name__ == "__main__":
    exp = 0 if len(sys.argv)<2 else int(sys.argv[1])
    main(exp)
