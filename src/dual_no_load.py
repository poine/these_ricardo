#-*- coding: utf-8 -*-

#
# Dual tethered planar VTOL
# See: https://poine.github.io/these_ricardo/planar_dual_no_load
#

import numpy as np, scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import misc_utils as mu


class Param:
    def __init__(self):
        self.m1, self.m2 = 1., 1.     # masses
        self.J1, self.J2 = 0.01, 0.01 # inertia
        self.d1, self.d2 = 0.2, 0.2   # half widths
        self.l = 1.                   # tether length
        self.g = 9.81                 # gravity
        self.compute_aux()

    def compute_aux(self):
        self.dovJ1, self.dovJ2 = self.d1/self.J1, self.d2/self.J2
        self.mt = self.m1+self.m2
        self.m1r = self.m1/self.mt
        self.m2r = self.m2/self.mt
        self.lm1r = self.l*self.m1r
        self.lm2r = self.l*self.m2r

class PVT:
    s_x, s_z, s_th, s_ph, s_th2, s_xd, s_zd, s_thd, s_phd, s_th2d, s_size  = range(11) # state vector components
    slice_pos, slice_kin, slice_dyn = slice(2), slice(5), slice(5, 11)
    i_fl1, i_fr1, i_fl2, i_fr2, i_size =  range(5)                                     # input vector components

    def __init__(self, P=None):
        self.P = P if P is not None else Param()
        self.Ue = self.P.g*np.array([self.P.m1/2, self.P.m1/2, self.P.m2/2, self.P.m2/2])

    def dyn(self, X, t, U):
        ut1, ut2 = (U[PVT.i_fl1]+U[PVT.i_fr1])/self.P.mt, (U[PVT.i_fl2]+U[PVT.i_fr2])/self.P.mt
        ud1, ud2 = (-U[PVT.i_fl1]+U[PVT.i_fr1])*self.P.dovJ1, (-U[PVT.i_fl2]+U[PVT.i_fr2])*self.P.dovJ2
        cph, sph = np.cos(X[PVT.s_ph]), np.sin(X[PVT.s_ph])
        ct1, st1 = np.cos(X[PVT.s_th]), np.sin(X[PVT.s_th])
        ct2, st2 = np.cos(X[PVT.s_th2]), np.sin(X[PVT.s_th2])

        if 0:
            A = np.array([[ 1,     0, -self.P.lm2r*sph],
                          [ 0,     1,  self.P.lm2r*cph],
                          [-sph, cph,  self.P.l]])
            invA = np.linalg.inv(A)

        
        
            xd, zd, phd = X[PVT.s_xd], X[PVT.s_zd], X[PVT.s_phd]
            ut2b = (U[PVT.i_fl1]+U[PVT.i_fr1])/self.P.m2
            B = np.array([[self.P.lm2r*phd**2*cph         -ut1*st1-ut2*st2],
                          [self.P.lm2r*phd**2*sph-self.P.g+ut1*ct1+ut2*ct2],
                          [-self.P.g*cph+ut2b*np.cos(X[PVT.s_ph]+X[PVT.s_th2])]])
            xdd, zdd, phdd =  (invA@B).flatten()
        else:
            phd = X[PVT.s_phd]
            phdd =  self.P.m1/self.P.l*ut2*np.cos(X[PVT.s_ph]-X[PVT.s_th2])-self.P.m2/self.P.l*ut1*np.cos(X[PVT.s_ph]-X[PVT.s_th])
            xdd = -ut1*st1 -ut2*st2+self.P.lm2r*(phdd*sph+phd**2*cph)
            zdd = ut1*ct1 +ut2*ct2 -self.P.g -self.P.lm2r*(phdd*cph-phd**2*sph)
            
        Xd = np.zeros(PVT.s_size)
        Xd[PVT.slice_kin] = X[PVT.slice_dyn]
        #breakpoint()
        Xd[PVT.s_xd:PVT.s_zd+1] = xdd, zdd
        Xd[PVT.s_phd] = phdd
        Xd[PVT.s_thd] = ud1
        Xd[PVT.s_th2d] = ud2
        return Xd

    def disc_dyn(self, Xk, Uk, dt):
        _unused, Xkp1 = scipy.integrate.odeint(self.dyn, Xk, [0, dt], args=(Uk,))
        return Xkp1

    def jac(self, Xe=None):
        Xe = Xe if Xe is not None else np.zeros(PVT.s_size)
        return mu.num_jacobian(Xe, self.Ue, self.dyn)

    def master_state(X): # x, z, theta, xd, zd, thd
        return np.hstack((X[:PVT.s_th+1], X[PVT.s_xd:PVT.s_thd+1]))
    def slave_state(X):  # phi, theta2, phid, theta2d
        return np.hstack((X[PVT.s_ph:PVT.s_th2+1], X[PVT.s_phd:PVT.s_th2d+1]))

    

def plot_trajectory(time, X, U, figure=None, axes=None):
    figure = plt.figure(tight_layout=True, figsize=[8., 6.]) if figure is None else figure
    nr, nc = (3,2) 
    axes = figure.subplots(nr, nc, sharex=True,) if axes is None else axes
    plots = [("$x$",             "m",     X[:,PVT.s_x]),
             ("$z$",             "m",     X[:,PVT.s_z]),
             ("$\\theta1$",      "deg",   np.rad2deg(X[:, PVT.s_th])),
             ("$\\theta2$",      "deg",   np.rad2deg(X[:, PVT.s_th2])),
             ("$\\phi$",         "deg",   np.rad2deg(X[:, PVT.s_ph]))]
    if U is not None:
        plots += [("$U$", "N", U)]
    for ax, (_t, _u, _d) in zip(axes.flatten(), plots):
        ax.plot(time, _d);  ax.set_title(_t); ax.grid(True); ax.yaxis.set_label_text(_u)
    axes[-1,0].xaxis.set_label_text('time in s'); axes[-1,1].xaxis.set_label_text('time in s')

    return figure, axes



#
# Animations
#
def extrema(arr, d=0): return np.min(arr)-d, np.max(arr)+d
def compute_extends(X, P):
    m_pos = X[:,PVT.slice_pos]
    s_pos = m_pos + P.P.l*np.array([np.cos(X[:,PVT.s_ph]), np.sin(X[:,PVT.s_ph])]).T
    (x_min, x_max), (y_min, y_max) = [extrema(np.hstack((m_pos[:,i], s_pos[:,i])), 1.2*P.P.d1) for i in range(2)]
    return (x_min, x_max), (y_min, y_max)
    
def animate(time, X, U, P, Xsp=None, figure=None, ax=None):
    (_xmin, _xmax), (_ymin, _ymax) = compute_extends(X, P)
    #_xmin, _xmax, _ymin, _ymax = -0.75, 1.75, -0.5, 0.5
    fig = figure or plt.figure(figsize=(10., 4.))
    if ax is None:
        ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(_xmin, _xmax),
                             ylim=(_ymin, _ymax), facecolor=(0.5, 0.9, 0.9))
    _time_template = 'time = {:0.1f}s'
    _time_text = ax.text(0.025, 0.92, 'N/A', transform=ax.transAxes)
    _line_body1, = ax.plot([], [], 'o-', lw=3, color='r', zorder=1)
    _line_body2, = ax.plot([], [], 'o-', lw=3, color='g', zorder=1)
    _line_tether, = ax.plot([], [], 'o-', lw=1, color='b', zorder=2)
    _line_sp, = ax.plot([], [], 'o--', lw=1, color='k', alpha=0.5, zorder=2)
    _lines = (_line_body1, _line_body2, _line_tether, _line_sp)
    
    def draw_quad(line, pos, theta, d):
        n = [d*np.cos(theta), d*np.sin(theta)]
        p1, p2 = pos + n, pos-n
        line.set_data([p1[0], p2[0]], [p1[1], p2[1]])
    
        
    def init():
        for _l in _lines: _l.set_data([], [])
        return (_time_text,) +  _lines

    def animate(i):
        p0 = X[i, P.slice_pos]
        n = np.array([np.cos(X[i, P.s_ph]), np.sin(X[i, P.s_ph])]) # tether aligned unit vector
        p1 = p0 + P.P.l*n
        _line_tether.set_data([p0[0], p1[0]], [p0[1], p1[1]])
        draw_quad(_line_body1, p0, X[i, P.s_th], P.P.d1 )
        draw_quad(_line_body2, p1, X[i, P.s_th2], P.P.d2 )
        if Xsp is not None:
            p0, phi0 = Xsp[i, P.slice_pos], Xsp[i, P.s_ph]
            n = np.array([np.cos(phi0), np.sin(phi0)]) # setpoint tether aligned unit vector
            p2 = p0 + P.P.l*n
            _line_sp.set_data([p0[0], p2[0]], [p0[1], p2[1]])
        _time_text.set_text(_time_template.format(i * dt))
        return (_time_text,) +  _lines

    dt = time[1]-time[0]; dt_mili = dt*1000
    anim = animation.FuncAnimation(fig, animate, np.arange(1, len(time)),
                                   interval=dt_mili, blit=True, init_func=init, repeat_delay=200)
    
    return anim
