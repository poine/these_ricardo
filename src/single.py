#!/usr/bin/env python
#-*- coding: utf-8 -*-

'''
  Dynamic model of a planar VTOL
  see: https://poine.github.io/these_ricardo/planar_single
'''

import numpy as np, matplotlib.pyplot as plt
import scipy.integrate

'''  Parameters of the dynamic model '''
class Param:
    def __init__(self, sat=None):
        self.m     = 1.      # mass in kg
        self.l     = 0.2     # arm length in m
        self.J     = 0.01    # inertia
        self.sat   = sat     # max thrust of motor in N
        self.g     = 9.81

''' State vector components '''
sv_x, sv_z, sv_th, sv_xd, sv_zd, sv_thd, sv_size = np.arange(0,7)
sv_slice_pos, sv_slice_vel = slice(2), slice(3,5)
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

def disc_dyn(Xk, Uk, P, dt):
    _unused, Xkp1 = scipy.integrate.odeint(dyn, Xk, [0, dt], args=(Uk,P))
    return Xkp1


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

#
# Animations
#
import matplotlib.animation as animation
def extrema(arr, d=0): return np.min(arr)-d, np.max(arr)+d
def compute_extends(X, P):
    m_pos = X[:,sv_slice_pos]
    (x_min, x_max), (y_min, y_max) = [extrema(m_pos[:,i], 1.5*P.l) for i in range(2)]
    return (x_min, x_max), (y_min, y_max)
def animate(time, X, U, Yc, P, title=None, _drawings=False, _imgs=True, figure=None, ax=None):
    (_xmin, _xmax), (_ymin, _ymax) = compute_extends(X, P)
    #_xmin, _xmax, _ymin, _ymax = -5, 5, -2, 2
    fig = figure or plt.figure(figsize=(10., 4.))
    if ax is None:
        ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(_xmin, _xmax),
                             ylim=(_ymin, _ymax), facecolor=(0.5, 0.9, 0.9))
    time_template = 'time = {:0.1f}s'
    time_text = ax.text(0.025, 0.92, 'Hello', transform=ax.transAxes)
    _line_body, = ax.plot([], [], 'o-', lw=3, color='r', zorder=1)
    _line_setpoint, = ax.plot([], [], 'o-', lw=1, color='b', zorder=1)

    #breakpoint()
    def init():
        #print('in init')
        _line_body.set_data([], [])
        return time_text, _line_body

    def animate(i):
        #print(f'in animate {i}')
        p0 = X[i, :sv_th]
        d = np.array([np.cos(X[i, sv_th]), np.sin(X[i, sv_th])])
        p1, p2 = p0+P.l*d, p0-P.l*d
        _line_body.set_data([p1[0], p0[0], p2[0]], [p1[1], p0[1], p2[1]])
        time_text.set_text(time_template.format(i * dt))
        p1 = p2 = Yc[i,0]
        _line_setpoint.set_data([p1[0], p2[0]], [p1[1], p2[1]])
        return time_text, _line_body, _line_setpoint
     
    dt = time[1]-time[0]
    dt_mili = dt*1000#25
    print(f'steps {len(time)} , {dt}, {dt_mili}')
    anim = animation.FuncAnimation(fig, animate, np.arange(1, len(time)),
                                   interval=dt_mili, blit=True, init_func=init, repeat_delay=200)

    return anim

