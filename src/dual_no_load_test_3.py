#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#
#
# Dual tethered, master-slave feedback regulation
#
#

import logging, sys
import numpy as np, scipy.integrate, matplotlib.pyplot as plt
import control.matlab

import dual_no_load as PVT, misc_utils as mu

LOG = logging.getLogger(__name__)

W = np.array([[1, 1],[-1, 1]])
invW = 0.5*np.array([[1, -1],[1, 1]]) # input variable change
def ms_feedback(X, Xsp, Ue, P):
    U = np.array(Ue)

    # Master
    omega_t, xi_t, omega_x, xi_x, omega_z, xi_z =6,0.9, 2,0.9, 2,0.9
    ox2, oz2, ot2 = omega_x**2, omega_z**2, omega_t**2
    txiomx, txiomz, txiomt = 2.*xi_x*omega_x, 2.*xi_z*omega_z, 2.*xi_t*omega_t  

    Kol = np.array([[0, oz2*P.P.m1, 0, txiomz*P.P.m1],
                    [-ox2/P.P.g, 0, -txiomx/P.P.g, 0]])
    Hol = np.array([[0, oz2*P.P.m1],[-ox2/P.P.g, 0]])
    Kil, Hil =  P.P.J1/P.P.d1*np.array([[ot2, txiomt]]), P.P.J1/P.P.d1*np.array([[ot2]])

    X_master = PVT.PVT.master_state(X) # x, z, th, xd, zd, thd
    # x, z, xd, zd
    X_master_ol = np.array([X_master[0], X_master[1], X_master[3], X_master[4]])
    Yc = Xsp[P.slice_pos]
    dut, dth = -np.dot(Kol, X_master_ol) + np.dot(Hol, Yc)
    # th, thd
    X_master_il =  np.array([X_master[2], X_master[5]])
    dud = -np.dot(Kil, X_master_il) + np.dot(Hil, dth)
    U[:2] += np.dot(invW, [dut, dud[0,0]])
    
    # Slave
    Kol = np.array([[0, oz2*P.P.m2, 0, txiomz*P.P.m2],
                    [-ox2/P.P.g, 0, -txiomx/P.P.g, 0]])
    Hol = np.array([[0, oz2*P.P.m2],[-ox2/P.P.g, 0]])
    Kil, Hil =  P.P.J2/P.P.d2*np.array([[ot2, txiomt]]), P.P.J2/P.P.d2*np.array([[ot2]])
    X_slave, X_slave_sp = P.slave_state2(X), P.slave_state2(Xsp)
    # x, z, xd, zd
    X_slave_ol = np.array([X_slave[0], X_slave[1], X_slave[3], X_slave[4]])
    Yc = X_slave_sp[:2]
    dut, dth = -np.dot(Kol, X_slave_ol) + np.dot(Hol, Yc)
    # th, thd
    X_slave_il =  np.array([X_slave[2], X_slave[5]])
    dud = -np.dot(Kil, X_slave_il) + np.dot(Hil, dth)
    U[2:] += np.dot(invW, [dut, dud[0,0]])

    return U

def stepx(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_x] = mu.step(t, 0.5)
    return X
def cst(t):
    X = np.zeros(PVT.PVT.s_size)
    X[PVT.PVT.s_x] = -0.25
    return X

def sim_feedback(sp, dt=0.01, save=None):
    P = PVT.PVT()
    X0, Ue = np.zeros(P.s_size), P.Ue

    X0 =np.array([0.1, 0.2, np.deg2rad(1),     np.deg2rad(1), np.deg2rad(2),
                  0.01, 0.02, np.deg2rad(0.1), np.deg2rad(3), np.deg2rad(4)])
    X0 =np.zeros(10)
    X0[PVT.PVT.s_x] = -0.1
    X0[PVT.PVT.s_z] = -0.1
    X0[PVT.PVT.s_zd] = 0.05
    X0[PVT.PVT.s_ph] = np.deg2rad(10)
    X0[PVT.PVT.s_phd] = -np.deg2rad(1)
    time = np.arange(0., 8., dt)
    X = np.zeros((len(time), P.s_size))
    Xsp = np.zeros((len(time), P.s_size))
    U = np.zeros((len(time), P.i_size))
    X[0] = X0
    for i in range(0, len(time)):
        Xsp[i] = sp(time[i])
        U[i] =  ms_feedback(X[i], Xsp[i], Ue, P)
        if i < len(time)-1: X[i+1] = P.disc_dyn(X[i], U[i], dt)
    PVT.plot_trajectory(time, X, U)
    if save is not None:
        plt.savefig(mu.PLOT_DIR+f'/dual_no_load__ms_ctl_{save}_chrono.png') 
    #anim = PVT.animate(time, X, U, P, Xsp)
    anim = PVT.Animation(time, X, U, P, Xsp).generate()

    if save is not None:
        mu.save_anim(mu.PLOT_DIR+f'/dual_no_load__ms_ctl_{save}.apng', anim)
    plt.show()
    
def main(exp, save):
    logging.basicConfig(level=logging.INFO)
    #sim_feedback(cst)
    sim_feedback(stepx, save='0')
    
if __name__ == "__main__":
    exp = 0 if len(sys.argv)<2 else int(sys.argv[1])
    save  = '-s' in sys.argv
    main(exp, save)
