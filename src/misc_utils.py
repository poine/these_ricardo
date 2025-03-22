import math, numpy as np


def step(t, a=1., p=8., dt=0.): return a if math.fmod(t+dt, p) > p/2 else -a

"""
Compute naïve numerical jacobian 
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




#  Linear reference models
#

# class LinRef:
#     ''' Linear Reference Model (with first order integration)'''
#     def __init__(self, K):
#         '''K: coefficients of the caracteristic polynomial, in ascending powers order,
#               highest order ommited (normalized to -1)'''
#         self.K = K; self.order = len(K)
#         self.X = np.zeros(self.order+1)

#     def run(self, dt, sp):
#         self.X[:self.order] += self.X[1:self.order+1]*dt
#         e =  np.array(self.X[:self.order]); e[0] -= sp
#         self.X[self.order] = np.sum(e*self.K)
#         return self.X

#     def poles(self):
#         return np.roots(np.insert(np.array(self.K[::-1]), 0, -1))

#     def reset(self, X0=None):
#         if X0 is None: X0 = np.zeros(self.order+1)
#         self.X = X0


# class FirstOrdLinRef(LinRef):
#     def __init__(self, tau):
#         LinRef.__init__(self, [-1/tau])

# class SecOrdLinRef(LinRef):
#     def __init__(self, omega, xi):
#         LinRef.__init__(self, [-omega**2, -2*xi*omega])



'''
  Linear reference models

 Adaptive Control with a Nested Saturation Reference Model
 https://pdfs.semanticscholar.org/8b5c/2718be69a651dc3233a934ac05f5edda7ffd.pdf
'''

class LinRef:
    ''' Nested Saturation Linear Reference Model (with first order integration)'''
    def __init__(self, K, sats=None):
        '''
        K: coefficients of the caracteristic polynomial, in ascending powers order,
              highest order ommited (normalized to -1)
        sats: saturations for each order in ascending order
        '''
        self.K = K; self.order = len(K)
        self.sats = sats
        if self.sats is not None:
            self.M = np.array(sats)
            self.M[0:-1] *= K[1:]
            for i in range(0, len(self.M)-1):
                self.M[len(self.M)-2-i] /= np.prod(self.M[len(self.M)-1-i:])
            self.CM = np.cumprod(self.M[::-1])[::-1]
            #print('M', self.M, 'CM', self.CM)
        self.X = np.zeros(self.order+1)

    def run(self, dt, sp):
        self.X[:self.order] += self.X[1:self.order+1]*dt
        e =  np.array(self.X[:self.order]); e[0] -= sp
        if self.sats is None:
            self.X[self.order] = np.sum(e*self.K)
        else:
            self.X[self.order] = 0
            for i in range(0, self.order):
                self.X[self.order] = self.M[i]*np.clip(self.K[i]/self.CM[i]*e[i] + self.X[self.order], -1., 1.)
        return self.X

    def poles(self):
        return np.roots(np.insert(np.array(self.K[::-1]), 0, -1))

    def reset(self, X0=None):
        if X0 is None: X0 = np.zeros(self.order+1)
        self.X = X0

class FirstOrdLinRef(LinRef):
    def __init__(self, tau):
        LinRef.__init__(self, [-1/tau])

class SecOrdLinRef(LinRef):
    def __init__(self, omega, xi, sats=None):
        LinRef.__init__(self, [-omega**2, -2*xi*omega], sats)

class ThirdOrderLinRef(LinRef):
    def __init__(self, omega, xi, tau, sats=None):
        LinRef.__init__(self, [-omega**2/tau, -1/tau-omega**2, -2*xi*omega-1/tau], sats)

def lamdba_to_omxi(l1, l2): return
def omxi_to_lambda(om, xi): re, im = -xi*om, np.sqrt(1-xi**2)*om; return complex(re,im), complex(re, -im)

# plotting helpers
PLOT_DIR='/home/poine/work/these_cogneti/these_ricardo/docs/plots'

def decorate(ax, title=None, xlab=None, ylab=None, legend=None, xlim=None, ylim=None, min_yspan=None):
    ax.xaxis.grid(color='k', linestyle='-', linewidth=0.2)
    ax.yaxis.grid(color='k', linestyle='-', linewidth=0.2)
    if xlab: ax.xaxis.set_label_text(xlab)
    if ylab: ax.yaxis.set_label_text(ylab)
    if title: ax.set_title(title, {'fontsize': 20 })

import time
import matplotlib.animation as animation
def save_anim(filename, an, dt=None):
    print('encoding animation video, please wait, it will take a while')
    _start = time.time()
    #fps = 1./dt/4; print(f'dt {dt} fps {fps}')# FIXME
    an.save(filename, writer=animation.PillowWriter(fps=90.)) # 90 85 80 75 is ok 91. 95. fails
    _end = time.time()
    print(f'video encoded, saved to {filename}, Bye (took {_end-_start:.1f} s)')
