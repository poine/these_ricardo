#!/usr/bin/env python
#-*- coding: utf-8 -*-

'''
  Dynamic model of a planar VTOL
  see: https://poine.github.io/ann_elucubrations/pvtol.html
'''

import numpy as np, matplotlib.pyplot as plt

'''  Parameters of the dynamic model '''
class Param:
    def __init__(self, sat=None):
        self.m     = 0.5     # mass in kg
        self.l     = 0.2     # length in m
        self.J     = 0.01    # inertia
        self.sat   = sat     # max thrust of motor in N
        self.g     = 9.81

''' State vector components '''
sv_x, sv_z, sv_th, sv_xd, sv_zd, sv_thd, sv_size = np.arange(0,7)
''' Input vector components '''
iv_f1, iv_f2, iv_size = np.arange(0,3)

def dyn(X, t, U, P):
    """
    State Space Representation: Xdot = f_param(t, X, U)
    """
    Xd = np.zeros(sv_size)
    Xd[sv_x:sv_th+1] = X[sv_xd:sv_thd+1]

    Usat = np.clip(U, 0, P.sat) if P.sat is not None else U
    ut = (Usat[iv_f1] + Usat[iv_f2])
    ud = (Usat[iv_f2] - Usat[iv_f1]) 

    sth, cth = np.sin(X[sv_th]), np.cos(X[sv_th])
    Xd[sv_xd]  = -sth/P.m*ut
    Xd[sv_zd]  =  cth/P.m*ut-P.g
    Xd[sv_thd] = P.l/P.J*ud

    return Xd    

def jacobian(Xe, Ue, P):
    """
    Partial derivatives of the dynamics:
    A = d/dX(dyn)(Xe, Ue), B = d/dU(dyn)(Xe, Ue)
    """   
    st, ct = np.sin(Xe[sv_th]), np.cos(Xe[sv_th])
    one_over_m = 1/P.m 
    A = np.array([[0, 0, 0, 1, 0, 0],
                  [0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 1],
                  [0, 0, -one_over_m*ct*(Ue[iv_f1]+Ue[iv_f2]), 0, 0, 0],
                  [0, 0, -one_over_m*st*(Ue[iv_f1]+Ue[iv_f2]),  0, 0, 0],
                  [0, 0, 0, 0, 0, 0]])
    l_over_J = P.l/P.J
    B = np.array([[0, 0],
                  [0, 0],
                  [0, 0],
                  [-one_over_m*st, -one_over_m*st],
                  [ one_over_m*ct,  one_over_m*ct],
                  [-l_over_J, l_over_J]])
    return A, B   



def trim(P, args={}):
    """
    Return equilibrium state Xe and input Ue such as
    dyn(Xe, Ue) = 0
    """
    Xe = np.zeros(sv_size)
    Xe[sv_x] = args.get('x', 0.)
    Xe[sv_z] = args.get('z', 0.)
    Ue =  0.5*P.m*9.81*np.ones(iv_size)
    return Xe, Ue



def plot_trajectory(time, X, U=None, figure=None, axes=None, window_title="Trajectory",
                    legend=None, filename=None):
    figure = plt.figure(tight_layout=True, figsize=[8., 6.]) if figure is None else figure
    nr, nc = (3,2) if U is None else (4,2)
    axes = figure.subplots(nr, nc, sharex=True,) if axes is None else axes
    plots = [("$x$",             "m",     X[:,sv_x]),
             ("$z$",             "m",     X[:,sv_z]),
             ("$\\dot{x}$",       "m/s",   X[:, sv_xd]),
             ("$\\dot{z}$",       "m/s",   X[:, sv_zd]),
             ("$\\theta$",       "deg",   np.rad2deg(X[:, sv_th])),
             ("$\\dot{\\theta}$", "deg/s", np.rad2deg(X[:, sv_thd]))]
    if U is not None:
        plots += [("$U$", "N", U)]
    for ax, (_t, _u, _d) in zip(axes.flatten(), plots):
        ax.plot(time, _d);  ax.set_title(_t); ax.grid(True); ax.yaxis.set_label_text(_u)
    axes[-1,0].xaxis.set_label_text('time in s'); axes[-1,1].xaxis.set_label_text('time in s')
    figure.canvas.manager.set_window_title(window_title)
    if filename is not None: plt.savefig(filename)
    return figure, axes



"""
Compute numerical jacobian 
"""
def num_jacobian(Xe, Ue, dyn):
    s_size = len(Xe)
    i_size = len(Ue)
    epsilonX = (0.1*np.ones(s_size)).tolist()
    dX = np.diag(epsilonX)
    A = np.zeros((s_size, s_size))
    for i in range(0, s_size):
        dx = dX[i,:]
        delta_f = dyn(Xe+dx/2, 0, Ue) - dyn(Xe-dx/2, 0, Ue)
        delta_f = delta_f / dx[i]
        A[:,i] = delta_f

    epsilonU = (0.1*np.ones(i_size)).tolist()
    dU = np.diag(epsilonU)
    B = np.zeros((s_size,i_size))
    for i in range(0, i_size):
        du = dU[i,:]
        delta_f = dyn(Xe, 0, Ue+du/2) - dyn(Xe, 0, Ue-du/2)
        delta_f = delta_f / du[i]
        B[:,i] = delta_f

    return A,B


