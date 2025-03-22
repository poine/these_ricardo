#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import sys
import numpy as np, matplotlib.pyplot as plt
# import matplotlib.animation as animation
# import matplotlib.image, matplotlib.offsetbox, matplotlib.transforms
# import time

import single as dyn
import single_test1
import misc_utils as mu
#
# Simple state feedback regulation
#


import control
def main(save=False):
    P = dyn.Param()
    Xe, Ue = dyn.trim(P)
    A, B = dyn.jacobian(Xe, Ue, P)
    Q, R = [5, 10, 1, 1, 1, 1], [1, 1]
    (K, X, E) = control.lqr(A, B, np.diag(Q), np.diag(R))
    K = np.array(K) # control.lqr returns a np.matrix object, we use np.array
    cl_dyn = A - np.dot(B,K)
    print('cl poles: {}'.format(np.linalg.eig(cl_dyn)[0]))
    C = np.array([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0]])
    H = single_test1.get_precommand(A, B, C, K)
    X0 = [1, 0, 0, 0, 0, 0]
    time, X, U = single_test1.sim_state_feedback(P, X0, Ue, K, H, single_test1.step_x, tf=8.)
    Yc = np.array([single_test1.step_x(_t) for _t in time]).reshape((len(time), 1, 2))
    #dyn.plot_trajectory(time, X, U, window_title="LQR")
    anim = dyn.animate(time, X, U, Yc, P)
    if save:
        mu.save_anim(mu.PLOT_DIR+'/single_step_x.apng', anim, 6/25.)
    plt.show()
  
if __name__ == "__main__":
    main(save='-s' in sys.argv)

