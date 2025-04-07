#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import sys
import numpy as np, matplotlib.pyplot as plt
# import matplotlib.animation as animation
# import matplotlib.image, matplotlib.offsetbox, matplotlib.transforms
# import time

import single as dyn
import misc_utils as mu
#
# Simple state feedback regulation
#

def stepx(t):
    X = np.zeros(dyn.sv_size)
    X[dyn.sv_x] = mu.step(t, 0.5)
    return X

def stepz(t):
    X = np.zeros(dyn.sv_size)
    X[dyn.sv_z] = mu.step(t, 0.5)
    X[dyn.sv_x] = mu.step(t, 0.5)
    return X

def circle(t, r=1, om=2.):
    X = np.zeros(dyn.sv_size)
    a = om*t
    X[dyn.sv_x] = r*np.cos(a)
    X[dyn.sv_z] = r*np.sin(a)
    return X

def sim_with_feedback(P, Ue, K, X0, sp, tf=10., dt=0.01):
    time = np.arange(0., tf, dt)
    X, Xsp, U = [np.zeros((len(time), _s)) for _s in (dyn.sv_size, dyn.sv_size, dyn.iv_size)]
    X[0] = X0
    for i in range(0, len(time)):
        Xsp[i] = sp(time[i])
        U[i] = Ue - K@(X[i]-Xsp[i])
        if i < len(time)-1: X[i+1] = dyn.disc_dyn(X[i], U[i], P, dt)
    return time, X, Xsp, U
    
import control
def main(exp, save=False):
    P = dyn.Param()
    Xe, Ue = dyn.trim(P)
    A, B = dyn.jacobian(Xe, Ue, P)
    Q, R = [10, 10, 1, 1, 1, 1], [1, 1]
    (K, X, E) = control.lqr(A, B, np.diag(Q), np.diag(R))
    K = np.array(K)
    print('cl poles: {}'.format(np.linalg.eig(A - B@K)[0]))
    if exp==0: sp, sf = stepx, 'single_step_x'
    elif exp==1: sp, sf = stepz, 'single_step_z'
    elif exp==2: sp, sf =lambda t: circle(t, om=0.5), 'single_circle_regulation'
    else: sp, sf = circle, 'single_circle_regulation'
    X0 = [1, 0, 0, 0, 0, 0]
    time, X, Xsp, U = sim_with_feedback(P, Ue, K, X0, sp, tf=10., dt=0.01)
    dyn.plot_trajectory(time, X, U, Xsp, window_title="LQR")
    if save:
        plt.savefig(mu.PLOT_DIR+f'/{sf}_chrono.png')
    anim = dyn.animate(time, X, U, P, Xsp[:,:2])
    if save:
        mu.save_anim(mu.PLOT_DIR+f'/{sf}.apng', anim)
    plt.show()
  
if __name__ == "__main__":
    exp = 0 if len(sys.argv)<2 else int(sys.argv[1])
    save  = '-s' in sys.argv
    main(exp, save)

