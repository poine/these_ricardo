#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import os, pdb
import math, numpy as np, matplotlib.pyplot as plt
import scipy.integrate
import control

import single as dyn

W = np.array([[1, 1],[-1, 1]])
invW = 0.5*np.array([[1, -1],[1, 1]])
def null_input_1(t): return np.zeros(1)
def null_input_2(t): return np.zeros(2)
def step(t, a=1., p=8., dt=0.): return a if math.fmod(t+dt, p) > p/2 else -a
def step_input_1(t): return np.array([step(t, np.deg2rad(10), dt=2)])
def step_x(t, a=1): return np.array([step(t, a), 0])
def step_z(t): return np.array([0, step(t)])
def step_xz(t): return np.array([step(t), step(t, dt=2)])


def sim_open_loop(P, X0, U):
    time = np.arange(0., 10, 0.01)
    X = scipy.integrate.odeint(dyn.dyn, X0, time, args=(U, P ))
    U = U*np.ones((len(time), dyn.iv_size))
    return time, X, U

def sim_state_feedback(P, X0, Ue, K, H, Ycf, tf=10.):
    def dU(X, t): return -np.dot(K, X) + np.dot(H, Ycf(t))
    def cl_dyn(X, t): return dyn.dyn(X, t, Ue+dU(X, t), P)
    time = np.arange(0., tf, 0.01)
    X = scipy.integrate.odeint(cl_dyn, X0, time)
    U =  np.array([Ue+dU(Xi, ti) for Xi, ti in zip(X, time)])
    return time, X, U


FIG_DIR = '/home/poine/work/tm_cmd_lin_multi/web/plots/td_2/'
#FIG_DIR = None
def savefig(fn):
    if FIG_DIR is not None and os.path.exists(FIG_DIR): plt.savefig(FIG_DIR+'/'+fn)

def question_1():
    P = dyn.Param()
    Ue = 0.5*P.m*P.g*np.ones(dyn.iv_size)
    X0 = [0, 0, 0, 0, 0, 0]
    time, X, U = sim_open_loop(P, X0, Ue)
    dyn.plot_trajectory(time, X, U, window_title='open loop equilibre'); savefig('01_open_loop_1.png')
    X0 = [0, 0, 0.1, 0, 0, 0]
    time, X, U = sim_open_loop(P, X0, Ue)
    dyn.plot_trajectory(time, X, U, window_title='open loop perturbé'); savefig('01_open_loop_2.png')

def question_2(omega=6, xi=0.9):
    P = dyn.Param()
    l_ov_J = P.l/P.J
    Ktilde = np.array([[0, 0, 0, 0, 0, 0],[0, 0, omega**2/l_ov_J, 0, 0, 2*xi*omega/l_ov_J]])
    Htilde = [[0], [omega**2/l_ov_J]]
    K, H = np.dot(invW, Ktilde), np.dot(invW, Htilde)
    Ue = 0.5*P.m*P.g*np.ones(dyn.iv_size)
    X0 = [0, 0, 0.1, 0, 0, 0]
    time, X, U = sim_state_feedback(P, X0, Ue, K, H, null_input_1)
    dyn.plot_trajectory(time, X, U, window_title='att stab feedback'); savefig('02_att_ctl_1.png')
    time, X, U = sim_state_feedback(P, X0, Ue, K, H, step_input_1)
    dyn.plot_trajectory(time, X, U, window_title='att stab setpoint'); savefig('02_att_ctl_2.png')

def get_precommand(A, B, C, K):
    tmp1 = np.linalg.inv(A - np.dot(B, K))
    tmp2 = np.dot(np.dot(C, tmp1), B)
    H = -np.linalg.inv(tmp2)
    return H

def test_rejection_perturbation(P, Ue, K, H, txt=''):
    for pert_X0, ax in zip([[0, 0, 0.1, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0]], ['theta', 'x', 'z']):
        time, X, U = sim_state_feedback(P, pert_X0, Ue, K, H, null_input_2, tf=4.)
        dyn.plot_trajectory(time, X, U, window_title="rejection perturbation {} {}".format(ax, txt))
        savefig(f'{txt}_rejection_pert.png')

def test_suivi_consigne(P, Ue, K, H, txt='', X0=np.zeros(dyn.sv_size)):
    for inp , w in zip([step_x, step_z, step_xz], ['x', 'z', 'xz']):
        time, X, U = sim_state_feedback(P, X0, Ue, K, H, inp)
        dyn.plot_trajectory(time, X, U, window_title="suivi consigne {} {}".format(w, txt))
        savefig(f'{txt}_suivi_consigne.png')

def test_variation_consigne_x(P, Ue, K, H, txt='', X0=np.zeros(dyn.sv_size)):
    for ampl in [1, 1.5, 2, 2.5]:
        time, X, U = sim_state_feedback(P, X0, Ue, K, H, lambda t: step_x(t, ampl))
        dyn.plot_trajectory(time, X, U, window_title="consigne x variable ampl {} {}".format(ampl, txt))  
        savefig(f'{txt}_consigne_variable.png')

def question_3(poles=[-20, -2-2j, -2+2j, -2.5, -2.5, -1.]):
    P = dyn.Param()
    Xe, Ue = dyn.trim(P)
    A, B = dyn.jacobian(Xe, Ue, P)
    K = np.array(control.place(A, B, poles))
    cl_dyn = A - np.dot(B,K)
    print('cl poles: {}'.format(np.linalg.eig(cl_dyn)[0]))
    C = np.array([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0]])
    H = get_precommand(A, B, C, K)
    test_rejection_perturbation(P, Ue, K, H, '03_Placement')
    test_suivi_consigne(P, Ue, K, H, '03_Placement')
    test_variation_consigne_x(P, Ue, K, H, '03_Placement')
      
def question_4(Q=[10, 1, 1, 1, 1, 1], R=[1, 1]):
    P = dyn.Param()
    Xe, Ue = dyn.trim(P)
    A, B = dyn.jacobian(Xe, Ue, P)
    (K, X, E) = control.lqr(A, B, np.diag(Q), np.diag(R))
    K = np.array(K) # control.lqr returns a np.matrix object, we use np.array
    cl_dyn = A - np.dot(B,K)
    print('cl poles: {}'.format(np.linalg.eig(cl_dyn)[0]))
    C = np.array([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0]])
    H = get_precommand(A, B, C, K)
    test_rejection_perturbation(P, Ue, K, H, '04_LQR')
    test_suivi_consigne(P, Ue, K, H, '04_LQR')
    test_variation_consigne_x(P, Ue, K, H, '04_LQR')

def question_5(omega_t=12, xi_t=0.9, omega_x=2, xi_x=0.9, omega_z=2, xi_z=0.9):
    P = dyn.Param()
    Xe, Ue = dyn.trim(P)
    l_ov_J = P.l/P.J
    ox2, oz2, ot2 = omega_x**2, omega_z**2, omega_t**2
    Ktilde = np.array([[0, oz2*P.m, 0, 0, 2*xi_z*omega_z*P.m, 0],
                       [-ot2*ox2/P.g/l_ov_J, 0, ot2/l_ov_J, -ot2*2*xi_x*omega_x/P.g/l_ov_J, 0, 2*xi_t*omega_t/l_ov_J]])
    Htilde = np.array([[0, oz2*P.m],[-ot2*ox2/P.g/l_ov_J, 0]])
    K, H = np.dot(invW, Ktilde), np.dot(invW, Htilde)
    A, B = dyn.jacobian(Xe, Ue, P)
    v, w = np.linalg.eig(A-np.dot(B, K)); print(v)
    test_rejection_perturbation(P, Ue, K, H, '05_NLI_total')
    test_suivi_consigne(P, Ue, K, H, '05_NLI_total')
    test_variation_consigne_x(P, Ue, K, H, '05_NLI_total')


def two_loop_ctl(X, t, Ue, Ycf, P, theta_sat=None):
    omega_t, xi_t, omega_x, xi_x, omega_z, xi_z =12,0.9, 2,0.9, 2,0.9
    ox2, oz2, ot2 = omega_x**2, omega_z**2, omega_t**2
    txiomx, txiomz, txiomt = 2.*xi_x*omega_x, 2.*xi_z*omega_z, 2.*xi_t*omega_t  
    Kol = np.array([[0, oz2*P.m, 0, txiomz*P.m],
                    [-ox2/P.g, 0, -txiomx/P.g, 0]])
    Hol = np.array([[0, oz2*P.m],[-ox2/P.g, 0]])
    Kil, Hil =  P.J/P.l*np.array([[ot2, txiomt]]), P.J/P.l*np.array([[ot2]])
    
    Xol = np.array([X[dyn.sv_x], X[dyn.sv_z], X[dyn.sv_xd], X[dyn.sv_zd]])
    dut, dth = -np.dot(Kol, Xol) + np.dot(Hol, Ycf(t))
    if theta_sat != None: dth = np.clip(dth, -theta_sat, theta_sat)
    Xil =  np.array([X[dyn.sv_th], X[dyn.sv_thd]])
    dud = -np.dot(Kil, Xil) + np.dot(Hil, dth)
    U = Ue + np.dot(invW, [dut, dud])
    return U

def question_6():
    P = dyn.Param(sat=6.)
    Xe, Ue = dyn.trim(P)
    def sim_with_ctl(theta_sat, X0, Ycf):
        def cl_dyn(X, t): return dyn.dyn(X, t, two_loop_ctl(X, t, Ue, Ycf, P, theta_sat), P)
        time = np.arange(0., 4., 0.01)
        X = scipy.integrate.odeint(cl_dyn, X0, time)
        U =  np.array([two_loop_ctl(Xi, ti, Ue, Ycf, P, theta_sat) for Xi, ti in zip(X, time)])
        return time, X, U

    X0, Ycf, fig, axes = [4, 0, 0, 0, 0, 0], null_input_2, None, None
    for theta_sat in [None, np.deg2rad(20), np.deg2rad(30), np.deg2rad(45)]:
        time, X, U =  sim_with_ctl(theta_sat, X0, Ycf)
        fig, axes = dyn.plot_trajectory(time, X, U, fig, axes)
    savefig('06_two_loops_theta_sat.png')
    
if __name__ == "__main__":
    #question_1()
    #question_2()
    #question_3()
    question_4()
    #question_5()
    #question_6()
    plt.show()
