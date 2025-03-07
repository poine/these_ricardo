#-*- coding: utf-8 -*-

import numpy as np, scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#from . import misc
import misc_utils as mu

'''

Dynamic model for the PVTOL pole system.
See https://poine.github.io/ann_elucubrations/pvtol_pole.html for equations

'''


class Param:
    def __init__(self):
        self.M = 1    # mass of the quad
        self.m = 0.5  # mass of the pole
        self.L = 0.5  # half width of the quand
        self.l = 0.5  # dist to pole cg
        self.J = 0.01 # inertia of the quad
        self.j = 0.01 # inertia of the pole
        self.g = 9.81
        self.compute_aux()
        
    def compute_aux(self):
        self.mt = self.M + self.m
        self.oneovmt = 1./self.mt
        self.mr = self.m / self.mt
        self.lmr = self.l*self.mr
        self.LovJ = self.L/self.J
        self.lmovjpml2 = self.l*self.m/(self.j+self.m*self.l*self.l)

class PVTP:
    s_x, s_z, s_th, s_ph, s_xd, s_zd, s_thd, s_phd, s_size = range(9) # state vector components
    i_f1, i_f2, i_size =  range(3) # input vector components
    
    def __init__(self, P=None):
        self.P = P if P is not None else Param()
        self.Ue = self.P.mt*self.P.g/2*np.ones(2)
    
    def dyn(self, X, t, U):
        sth, cth = np.sin(X[PVTP.s_th]), np.cos(X[PVTP.s_th])
        sph, cph = np.sin(X[PVTP.s_ph]), np.cos(X[PVTP.s_ph])
        phdsq = X[PVTP.s_phd]*X[PVTP.s_phd] # phi_dot_squared
        ut =  U[PVTP.i_f1] + U[PVTP.i_f2]
        ud = -U[PVTP.i_f1] + U[PVTP.i_f2]

        thdd = self.P.LovJ*ud # quad angular accel
        
        a = -self.P.lmr*cph
        e = -self.P.lmr*sph*phdsq -self.P.oneovmt*sth*ut
        b = -self.P.lmr*sph
        f =  self.P.lmr*cph*phdsq - self.P.g + self.P.oneovmt*cth*ut
        c = -self.P.lmovjpml2*cph
        d = -self.P.lmovjpml2*sph
        h =  self.P.lmovjpml2*sph*self.P.g
        A = np.array([[1, 0, a],[0, 1, b],[c, d, 1]])
        Y = np.array([[e],[f],[h]])
        xdd, zdd, phdd = np.dot(np.linalg.inv(A), Y)
        #breakpoint()
        Xd = np.zeros(8)
        Xd[:4] = X[4:]
        Xd[4:] = np.array([xdd, zdd, [thdd], phdd]).flatten()
        return Xd

    def jac(self):
        Xe = np.zeros(8)
        return mu.num_jacobian(Xe, self.Ue, self.dyn)

    def disc_dyn(self, Xk, Uk, dt):
        _unused, Xkp1 = scipy.integrate.odeint(self.dyn, Xk, [0, dt], args=(Uk,))
        return Xkp1


def plot_trajectory(time, X, U, figure=None, axes=None):
    figure = plt.figure(tight_layout=True, figsize=[8., 6.]) if figure is None else figure
    nr, nc = (4,2) if U is None else (5,2)
    axes = figure.subplots(nr, nc, sharex=True,) if axes is None else axes
    plots = [("$x$",             "m",     X[:,PVTP.s_x]),
             ("$z$",             "m",     X[:,PVTP.s_z]),
             ("$\\dot{x}$",       "m/s",   X[:, PVTP.s_xd]),
             ("$\\dot{z}$",       "m/s",   X[:, PVTP.s_zd]),
             ("$\\theta$",       "deg",   np.rad2deg(X[:, PVTP.s_th])),
             ("$\\dot{\\theta}$", "deg/s", np.rad2deg(X[:, PVTP.s_thd])),
             ("$\\phi$",       "deg",   np.rad2deg(X[:, PVTP.s_ph])),
             ("$\\dot{\\phi}$", "deg/s", np.rad2deg(X[:, PVTP.s_phd]))]
    if U is not None:
        plots += [("$U$", "N", U)]
    for ax, (_t, _u, _d) in zip(axes.flatten(), plots):
        ax.plot(time, _d);  ax.set_title(_t); ax.grid(True); ax.yaxis.set_label_text(_u)
    axes[-1,0].xaxis.set_label_text('time in s'); axes[-1,1].xaxis.set_label_text('time in s')
    return figure, axes



def animate(time, X, U, P, figure=None, ax=None):
    _xmin, _xmax = -5, 5
    _ymin, _ymax = -2, 2
    fig = figure or plt.figure(figsize=(10., 5.))
    if ax is None:
        ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(_xmin, _xmax),
                             ylim=(_ymin, _ymax), facecolor=(0.5, 0.9, 0.9))
    time_template = 'time = {:0.1f}s'
    time_text = ax.text(0.025, 0.92, 'Hello', transform=ax.transAxes)
    _line_body, = ax.plot([], [], 'o-', lw=3, color='r', zorder=1)
    
    def init():
        #print('in init')
        _line_body.set_data([], [])
        return time_text, _line_body
    def animate(i):
        #print(f'in animate {i}')
        p0 = X[i, :P.s_th]
        d1 = np.array([np.cos(X[i, P.s_th]), np.sin(X[i, P.s_th])])
        p1, p2 = p0+P.P.L*d1, p0-P.P.L*d1
        d2 = np.array([np.sin(X[i, P.s_ph]), np.cos(X[i, P.s_ph])])
        p3 = p0+2*P.P.l*d2
        _line_body.set_data([p1[0], p0[0], p3[0], p0[0], p2[0]], [p1[1], p0[1], p3[1], p0[1], p2[1]])
        time_text.set_text(time_template.format(i * dt))
        return time_text, _line_body
     
    dt = time[1]-time[0]
    dt_mili = dt*1000#25
    print(f'steps {len(time)} , {dt}, {dt_mili}')
    anim = animation.FuncAnimation(fig, animate, np.arange(1, len(time)),
                                   interval=dt_mili, blit=True, init_func=init, repeat_delay=200)

    return anim
