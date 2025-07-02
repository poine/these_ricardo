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

def normMpiPi(_v): return np.arctan2(np.sin(_v), np.cos(_v))
def normMpiPi(_v):
    _v = np.fmod(_v, 2*np.pi)
    if _v > np.pi: _v -= np.pi
    return _v

def cst(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_ph] = np.pi/3
    return X

def stepx(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_x] = mu.step(t, 0.5)
    X[PVT.PVT.s_ph] = np.pi/12
    return X

def stepz(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_z] = mu.step(t, 0.4)
    return X

def stepxz(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_x] = mu.step(t, 0.5, dt=2.)
    X[PVT.PVT.s_z] = mu.step(t, 0.4)
    return X

def stepphi(t, phi0=np.deg2rad(90.), dphi=np.deg2rad(15)):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_ph] = phi0 + mu.step(t, dphi)
    return X

def rampphi(t, omega=2*np.pi/10.):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_ph] = omega*t#np.sin(0.5*t)#
    return X

def step_foo(t, dphi=np.deg2rad(45)):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_x] = mu.step(t, 0.5)
    X[PVT.PVT.s_z] = -0.5
    X[PVT.PVT.s_ph] = np.pi/2 + mu.step(t, dphi)
    return X

def feedback_control(X, Xsp, Ue, K):
    dX = X-Xsp
    dX[PVT.PVT.s_ph] = normMpiPi(dX[PVT.PVT.s_ph])
    dX[PVT.PVT.s_th] = normMpiPi(dX[PVT.PVT.s_th])
    dX[PVT.PVT.s_th2] = normMpiPi(dX[PVT.PVT.s_th2])
    return Ue - K@dX


def sim_feedback(sp, dt=0.01, savefile=None):
    P = PVT.PVT()
    Q = np.diag([3., 7.,10., 0.5, 2.,
                 0.25, 0.25, 0.25, 0.005, 0.01])
    R = np.diag([0.5, 0.5, 0.5, 0.5])
    def compute_gain(_X, _U, _v=False):
        _Xe = np.array(_X)
        _Xe[P.s_th] = _Xe[P.s_th2] = 0
        _Xe[P.slice_dyn] = 0
        A, B = P.jac(_Xe)
        (K, X, E) = control.matlab.lqr(A, B, Q, R)
        if _v:
            with np.printoptions(precision=2, linewidth=200):
                poles, vect_p = np.linalg.eig(A-B@K)
                LOG.info('gain:\n{}'.format(K))
                LOG.info('poles:\n{}'.format(poles))
        return K
    X0, Ue = np.zeros(P.s_size), P.Ue
    #K = compute_gain(X0, Ue)
    X0[P.s_th] = np.deg2rad(11.5)
    X0[P.s_th2] = np.deg2rad(-11.5)
    #X0[P.s_ph] = 2*np.pi/3
    time = np.arange(0., 12.7, dt)
    X, Xsp, U = [np.zeros((len(time), _s)) for _s in (P.s_size, P.s_size, P.i_size)]
    X[0] = X0
    for i in range(0, len(time)):
        Xsp[i] = sp(time[i])
        K = compute_gain(X[i], P.Ue) 
        U[i] =  feedback_control(X[i], Xsp[i], Ue, K)
        if i < len(time)-1: X[i+1] = P.disc_dyn(X[i], U[i], dt)
    PVT.plot_trajectory(time, X, U)
    anim = PVT.Animation(time, X, U, P, Xsp).generate()

    if savefile is not None:
        mu.save_anim(mu.PLOT_DIR+'/'+savefile, anim)
    plt.show()


def test_0(save): sim_feedback(cst)
def test_1(save): sim_feedback(stepx, savefile='dual_no_load__state_feedback_1.apng' if save else None)
def test_2(save): sim_feedback(stepz)
def test_3(save): sim_feedback(stepphi)
def test_4(save): sim_feedback(stepxz, savefile='dual_no_load__state_feedback_2.apng' if save  else None)
def test_5(save): sim_feedback(rampphi, savefile='dual_no_load__state_feedback_3.apng' if save else None)
def test_6(save): sim_feedback(step_foo, savefile='dual_no_load__state_feedback_4.apng' if save else None)

def main(exp, save):
    logging.basicConfig(level=logging.INFO)
    exps = [test_0, test_1, test_2, test_3, test_4, test_5, test_6]
    exps[exp](save)

if __name__ == "__main__":
    exp = 1 if len(sys.argv)<2 else int(sys.argv[1])
    save  = '-s' in sys.argv
    main(exp, save)
