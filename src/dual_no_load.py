#-*- coding: utf-8 -*-

#
# Dual tethered planar VTOL
# See: https://poine.github.io/these_ricardo/planar_dual_no_load
#

import numpy as np, scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import misc_utils as mu, anim_utils

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
    s_x, s_z, s_ph, s_th, s_th2, s_xd, s_zd, s_phd, s_thd, s_th2d, s_size  = range(11) # state vector components
    slice_pos, slice_vel, slice_kin, slice_dyn = slice(2), slice(5, 7), slice(5), slice(5, 11)
    i_fl1, i_fr1, i_fl2, i_fr2, i_size =  range(5)                                     # input vector components

    def __init__(self, P=None):
        self.P = P if P is not None else Param()
        self.Ue = self.P.g*np.array([self.P.m1/2, self.P.m1/2, self.P.m2/2, self.P.m2/2])

    def dyn_wrong(self, X, t, U):
        ut1, ut2 = (U[PVT.i_fl1]+U[PVT.i_fr1])/self.P.mt, (U[PVT.i_fl2]+U[PVT.i_fr2])/self.P.mt
        ud1, ud2 = (-U[PVT.i_fl1]+U[PVT.i_fr1])*self.P.dovJ1, (-U[PVT.i_fl2]+U[PVT.i_fr2])*self.P.dovJ2
        cph, sph = np.cos(X[PVT.s_ph]), np.sin(X[PVT.s_ph])
        ct1, st1 = np.cos(X[PVT.s_th]), np.sin(X[PVT.s_th])
        ct2, st2 = np.cos(X[PVT.s_th2]), np.sin(X[PVT.s_th2])

        phd = X[PVT.s_phd]
        phdd =  self.P.m1/self.P.l*ut2*np.cos(X[PVT.s_ph]-X[PVT.s_th2])-self.P.m2/self.P.l*ut1*np.cos(X[PVT.s_ph]-X[PVT.s_th])
        xdd = -ut1*st1 -ut2*st2+self.P.lm2r*(phdd*sph+phd**2*cph)
        zdd = ut1*ct1 +ut2*ct2 -self.P.g -self.P.lm2r*(phdd*cph-phd**2*sph)
            
        Xd = np.zeros(PVT.s_size)
        Xd[PVT.slice_kin] = X[PVT.slice_dyn]
        
        Xd[PVT.s_xd:PVT.s_zd+1] = xdd, zdd
        Xd[PVT.s_phd] = phdd
        Xd[PVT.s_thd] = ud1
        Xd[PVT.s_th2d] = ud2
        return Xd

    def dyn(self, X, t, U):
        m1, m2, mt, l, g = self.P.m1, self.P.m2, self.P.mt, self.P.l, self.P.g
        F1, F2 = U[0]+U[1], U[2]+U[3]
        ph, th1, th2, phd = X[PVT.s_ph], X[PVT.s_th], X[PVT.s_th2], X[PVT.s_phd]
        phmth1, phmth2 = ph-th1, ph-th2
        cphmth1, sphmth1 = np.cos(phmth1), np.sin(phmth1)
        cphmth2, sphmth2 = np.cos(phmth2), np.sin(phmth2)
        cph, sph = np.cos(X[PVT.s_ph]), np.sin(X[PVT.s_ph])
        cth1, sth1 = np.cos(X[PVT.s_th]), np.sin(X[PVT.s_th])
        cth2, sth2 = np.cos(X[PVT.s_th2]), np.sin(X[PVT.s_th2])
        #xdd = (-F1*sth1 -F2*sth2 +l*m2*cph*phd**2 + (F1*m2*cphmth1-F2*m1*cphmth2)*sph/m2)/mt
        #zdd = (F1*cth1 + F2*cth2  -g*mt +...TODO
        phdd = (F1*m2*cphmth1 - F2*m1*cphmth2)/(l*m2**2)
        # version using phdd
        xdd = (-F1*sth1 -F2*sth2 +l*m2*(sph*phdd + cph*phd**2))/mt
        zdd = (F1*cth1 + F2*cth2 -g*mt +l*m2*sph*phd**2 -l*m2*cph*phdd)/mt
        th1dd = (-U[0]+U[1])*self.P.dovJ1
        th2dd = (-U[2]+U[3])*self.P.dovJ2
        
        Xd = np.zeros(PVT.s_size)
        Xd[PVT.slice_kin] = X[PVT.slice_dyn]
        Xd[PVT.slice_dyn] = xdd, zdd, phdd, th1dd, th2dd
        return Xd
        
    
    def disc_dyn(self, Xk, Uk, dt):
        _unused, Xkp1 = scipy.integrate.odeint(self.dyn, Xk, [0, dt], args=(Uk,))
        return Xkp1

    def jac(self, Xe=None):
        Xe = Xe if Xe is not None else np.zeros(PVT.s_size)
        return mu.num_jacobian(Xe, self.Ue, self.dyn)

    # FIXME orders
    # def master_state(X): # x, z, theta, xd, zd, thd
    #     return np.hstack((X[:PVT.s_th+1], X[PVT.s_xd:PVT.s_thd+1]))
    # def slave_state(X):  # phi, theta2, phid, theta2d
    #     return np.hstack((X[PVT.s_ph:PVT.s_th2+1], X[PVT.s_phd:PVT.s_th2d+1]))

    # def slave_state2(self, X):
    #     #
    #     cph, sph, phd = np.cos(X[PVT.s_ph]), np.sin(X[PVT.s_ph]), X[PVT.s_phd]
    #     x2, z2 = X[PVT.slice_pos] + self.P.l*np.array([cph, sph])
    #     x2d, z2d = X[PVT.slice_vel] + self.P.l*phd*np.array([-sph, cph])                         
    #     return np.array([x2, z2, X[PVT.s_th2], x2d, z2d, X[PVT.s_th2d]])

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

class Animation(anim_utils.Animation):
    def setup_drawing(self):
        self._line_body1, = self.ax.plot([], [], 'o-', lw=3, color='r', zorder=1)
        self._line_body2, = self.ax.plot([], [], 'o-', lw=3, color='g', zorder=1)
        self._line_tether, = self.ax.plot([], [], 'o-', lw=1, color='b', zorder=2)
        self._line_sp, = self.ax.plot([], [], 'o--', lw=1, color='k', alpha=0.5, zorder=2)
        self._lines = (self._line_body1, self._line_body2, self._line_tether, self._line_sp)
        return self._lines
        
    def animate_drawing(self, i):
        p0 = self.X[i, PVT.slice_pos] # center of quad 1
        n = np.array([np.cos(self.X[i, PVT.s_ph]), np.sin(self.X[i, PVT.s_ph])]) # tether aligned unit vector
        p1 = p0 + self.P.P.l*n # center of quad 2
        self._line_tether.set_data([p0[0], p1[0]], [p0[1], p1[1]])
        self.draw_quad(self._line_body1, p0, self.X[i, PVT.s_th], self.P.P.d1 )
        self.draw_quad(self._line_body2, p1, self.X[i, PVT.s_th2], self.P.P.d2 )
        if self.Xsp is not None:
            p0, phi0 = self.Xsp[i, PVT.slice_pos], self.Xsp[i, PVT.s_ph]
            n = np.array([np.cos(phi0), np.sin(phi0)]) # setpoint tether aligned unit vector
            p2 = p0 + self.P.P.l*n
            self._line_sp.set_data([p0[0], p2[0]], [p0[1], p2[1]])
        return self._lines

    def compute_extends(self, X, P):
        m_pos = self.X[:,PVT.slice_pos]
        s_pos = m_pos + self.P.P.l*np.array([np.cos(X[:,PVT.s_ph]), np.sin(X[:,PVT.s_ph])]).T
        (x_min, x_max), (y_min, y_max) = [anim_utils.extrema(np.hstack((m_pos[:,i], s_pos[:,i])), 1.2*P.P.d1) for i in range(2)]
        return (x_min, x_max), (y_min, y_max)
         
