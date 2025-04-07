#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import sys, numpy as np, scipy.integrate
import matplotlib.pyplot as plt
import control.matlab

import single_point_load as dyn, misc_utils as mu

def sim_open_loop(P, X0, U):
    time = np.arange(0., 10, 0.01)
    X = scipy.integrate.odeint(P.dyn, X0, time, args=(U, ))
    U = U*np.ones((len(time), P.i_size))
    return time, X, U

def test_open_loop(save=False, sf='single_point_load__open_loop'):
    P = dyn.PVTP()
    X0 = np.zeros(P.s_size)
    X0[P.s_th] = np.deg2rad(0.5 )
    X0[P.s_ph] = np.deg2rad(-1.)
    time, X, U = sim_open_loop(P, X0, P.Ue)
    dyn.plot_trajectory(time, X, U)
    if save:
        plt.savefig(mu.PLOT_DIR+f'/{sf}_chrono.png')
        anim = dyn.animate(time, X, U, P)
    if save:
        mu.save_anim(mu.PLOT_DIR+f'/{sf}.apng', anim)
    plt.show()

    
if __name__ == "__main__":
    test_open_loop(save='-s' in sys.argv)
