#!/usr/bin/env python3
import sys
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt

'''
  Linear reference models

 Adaptive Control with a Nested Saturation Reference Model
 https://pdfs.semanticscholar.org/8b5c/2718be69a651dc3233a934ac05f5edda7ffd.pdf
'''

import misc_utils as mu

def plot_ref(_time, Xr, _id="", _sp=None, _sats=None, fig=None, axs=None):
    _, order = Xr.shape
    fig = plt.figure(tight_layout=True, figsize=[8., 2.*order]) if fig is None else fig
    axs = fig.subplots(order, 1, sharex=True,) if axs is None else axs
    labels = ['pos', 'vel', 'accel', 'jerk', 'snap', 'crackle', 'pop']
    for _i, (_a, _d, _l) in enumerate(zip(axs, Xr.T, labels)):
        _a.plot(_time, _d, label=_id)
        if _sats is not None and _i>0:
            _a.plot(_time, _sats[_i-1]*np.ones(len(_time)), 'k--')
            _a.plot(_time, -_sats[_i-1]*np.ones(len(_time)), 'k--')
        mu.decorate(_a, _l, xlab=None if _i<order-1 else 'time in s', legend=True)
    if _sp is not None: axs[0].plot(_time, _sp, label='setpoint')
    return fig, axs

def run_ref(_ref, _time, _sp):
    Xr = np.zeros((len(_time), _ref.order+1))
    for i in range(1, len(_time)):
        dt = _time[i] - _time[i-1]
        Xr[i] = _ref.run(dt, _sp[i])
    return Xr

def test_2nd_order():
    _time = np.arange(0, 10, 0.01)

    _sats = [6., 35.]  # vel, accel
    _ref1 = mu.SecOrdLinRef(omega=6, xi=0.9, sats=None)
    _ref2 = mu.SecOrdLinRef(omega=6, xi=0.9, sats=_sats)

    _sp = 4.*scipy.signal.square(_time*np.pi/3)
    Xr1 = run_ref(_ref1, _time, _sp)
    Xr2 = run_ref(_ref2, _time, _sp)

    fig, axs = plot_ref(_time, Xr1, _id="linear")
    plot_ref(_time, Xr2, _id="saturated", _sp=_sp, _sats=_ref2.sats, fig=fig, axs=axs)
    plt.show()

def test_5th_order():
    _time = np.arange(0, 12, 0.01) 

    poles = [complex(-3, 3), complex(-3, -3), complex(-4, 4), complex(-4, -4), -5]
    coefs = np.poly(poles)
    _ref1 = mu.LinRef(-np.flip(coefs)[:-1], sats=None) 
    _sats = [7., 15., 100., 500, 15000.]  # vel, accel, jerk, snap, crackle 
    _ref2 = mu.LinRef(-np.flip(coefs)[:-1], sats=_sats)

    _sp = 10.*scipy.signal.square(_time*np.pi/5)
    Xr1 = run_ref(_ref1, _time, _sp)
    Xr2 = run_ref(_ref2, _time, _sp)
    fig, axs = plot_ref(_time, Xr1, _id="linear")
    plot_ref(_time, Xr2, _id="saturated", _sp=_sp, _sats=_ref2.sats, fig=fig, axs=axs)
    plt.show()

def test_6th_order():
    _time = np.arange(0, 12, 0.01) 

    poles = [complex(-3, 3), complex(-3, -3), complex(-4, 4), complex(-4, -4), complex(-5, 5), complex(-5, -5)]
    coefs = np.poly(poles)
    _ref1 = mu.LinRef(-np.flip(coefs)[:-1], sats=None) 
    _sats = [7., 15., 100., 500, 15000., 100000]  # vel, accel, jerk, snap, crackle, pop 
    _ref2 = mu.LinRef(-np.flip(coefs)[:-1], sats=_sats)

    _sp = 10.*scipy.signal.square(_time*np.pi/5)
    Xr1 = run_ref(_ref1, _time, _sp)
    Xr2 = run_ref(_ref2, _time, _sp)
    fig, axs = plot_ref(_time, Xr1, _id="linear")
    plot_ref(_time, Xr2, _id="saturated", _sp=_sp, _sats=_ref2.sats, fig=fig, axs=axs)
    plt.show()
        
def main(args):
    #test_2nd_order()
    #test_5th_order()
    test_6th_order()

if __name__ == '__main__':
    main(sys.argv)
