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
        self.d     = 0.2     # arm length in m
        self.J     = 0.01    # inertia in kg⋅m2
        self.sat   = sat     # max thrust of motor in N
        self.g     = 9.81
        self.compute_aux()
        
    def compute_aux(self):
        self.one_over_m = 1./self.m
        self.d_over_J = self.d/self.J
        
''' State vector components '''
sv_x, sv_z, sv_th, sv_xd, sv_zd, sv_thd, sv_size = np.arange(0,7)
sv_slice_pos, sv_slice_vel = slice(2), slice(3,5)
sv_lice_kin, sv_slice_dyn = slice(3), slice(3,6)
''' Input vector components '''
iv_fl, iv_fr, iv_size = np.arange(0,3)

def dyn(X, t, U, P):
    """
    State Space Representation: Xdot = f_param(t, X, U)
    """
    Xd = np.zeros(sv_size)
    Xd[sv_x:sv_th+1] = X[sv_xd:sv_thd+1]

    Usat = np.clip(U, 0, P.sat) if P.sat is not None else U
    ut = ( Usat[iv_fl] + Usat[iv_fr])
    ud = (-Usat[iv_fl] + Usat[iv_fr]) 

    sth, cth = np.sin(X[sv_th]), np.cos(X[sv_th])
    Xd[sv_xd]  = -sth/P.m*ut
    Xd[sv_zd]  =  cth/P.m*ut-P.g
    Xd[sv_thd] = P.d/P.J*ud

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
    A = np.array([[0, 0, 0, 1, 0, 0],
                  [0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 1],
                  [0, 0, -P.one_over_m*ct*(Ue[iv_fl]+Ue[iv_fr]), 0, 0, 0],
                  [0, 0, -P.one_over_m*st*(Ue[iv_fl]+Ue[iv_fr]), 0, 0, 0],
                  [0, 0, 0, 0, 0, 0]])
    B = np.array([[0, 0],
                  [0, 0],
                  [0, 0],
                  [-P.one_over_m*st, -P.one_over_m*st],
                  [ P.one_over_m*ct,  P.one_over_m*ct],
                  [-P.d_over_J, P.d_over_J]])
    return A, B   



def trim(P, args={}):
    """
    Return equilibrium state Xe and input Ue such as
    dyn(Xe, Ue) = 0
    """
    Xe = np.zeros(sv_size)
    Xe[sv_x], Xe[sv_z] = args.get('x', 0.), args.get('z', 0.)
    Ue =  0.5*P.m*P.g*np.ones(iv_size)
    return Xe, Ue



def plot_trajectory(time, X, U=None, Yr=None, Ysp=None, Xr=None,
                    figure=None, axes=None, window_title="Trajectory",
                    legend=None, filename=None):
    figure = plt.figure(tight_layout=True, figsize=[8., 6.]) or figure
    nr, nc = (3,2) if U is None else (4,2)
    axes = figure.subplots(nr, nc, sharex=True,) if axes is None else axes
    plots = [("$x$",              "m",     X[:,sv_x]),
             ("$z$",              "m",     X[:,sv_z]),
             ("$\\dot{x}$",       "m/s",   X[:, sv_xd]),
             ("$\\dot{z}$",       "m/s",   X[:, sv_zd]),
             ("$\\theta$",        "deg",   np.rad2deg(X[:, sv_th])),
             ("$\\dot{\\theta}$", "deg/s", np.rad2deg(X[:, sv_thd]))]
    if U is not None:
        plots += [("$U$", "N", U)]
    for ax, (_t, _u, _d) in zip(axes.flatten(), plots):
        ax.plot(time, _d);  ax.set_title(_t); ax.grid(True); ax.yaxis.set_label_text(_u)
    if Xr is not None:
        axes[0,0].plot(time, Xr[:,0], label='reference') 
        axes[0,1].plot(time, Xr[:,1], label='reference') 
        axes[1,0].plot(time, Xr[:,3], label='reference') 
        axes[1,1].plot(time, Xr[:,4], label='reference') 
        axes[2,0].plot(time, np.rad2deg(Xr[:,2]), label='reference') 
        axes[2,1].plot(time, np.rad2deg(Xr[:,5]), label='reference') 
    if Ysp is not None:
        axes[0,0].plot(time, Ysp[:,0], label='setpoint') 
        axes[0,1].plot(time, Ysp[:,1], label='setpoint') 
    if Yr is not None:
        axes[0,0].plot(time, Yr[:,0], label='reference') 
        axes[0,1].plot(time, Yr[:,1], label='reference') 
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
    (x_min, x_max), (y_min, y_max) = [extrema(m_pos[:,i], 1.5*P.d) for i in range(2)]
    return (x_min, x_max), (y_min, y_max)
def subsample(l): # FIXME hardcoded 100Hz -> 25fps 
    return [_l[::4] if _l is not None else None for _l in l]

def animate(time, X, U, P, Ysp=None, Yr=None, title=None, _drawings=False, _imgs=True, figure=None, ax=None):
    time, X, U, Ysp, Yr =  subsample((time, X, U, Ysp, Yr))
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
    _line_reference, = ax.plot([], [], 'o-', lw=1, color='g', zorder=1)
    
    #breakpoint()
    def init():
        #print('in init')
        _line_body.set_data([], [])
        return time_text, _line_body

    def animate(i):
        #print(f'in animate {i}')
        p0 = X[i, :sv_th]
        d = np.array([np.cos(X[i, sv_th]), np.sin(X[i, sv_th])])
        p1, p2 = p0+P.d*d, p0-P.d*d
        _line_body.set_data([p1[0], p0[0], p2[0]], [p1[1], p0[1], p2[1]])
        time_text.set_text(time_template.format(i * dt))
        if Ysp is not None:
            p1 = p2 = Ysp[i]
            _line_setpoint.set_data([p1[0], p2[0]], [p1[1], p2[1]])
        if Yr is not None:
            p1 = p2 = Yr[i]
            _line_reference.set_data([p1[0], p2[0]], [p1[1], p2[1]])
        return time_text, _line_body, _line_setpoint, _line_reference
     
    dt = time[1]-time[0]
    dt_mili = dt*1000#25
    print(f'steps {len(time)} , {dt}, {dt_mili}')
    anim = animation.FuncAnimation(fig, animate, np.arange(1, len(time)),
                                   interval=dt_mili, blit=True, init_func=init, repeat_delay=200)

    return anim

