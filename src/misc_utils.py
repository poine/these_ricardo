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

