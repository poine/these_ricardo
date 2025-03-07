#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import numpy as np, matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.image, matplotlib.offsetbox, matplotlib.transforms
import time

import single as dyn
import single_test1

def animate(time, X, U, Yc, P, title=None, _drawings=False, _imgs=True, figure=None, ax=None):
    _xmin, _xmax = -5, 5
    _ymin, _ymax = -2, 2
    fig = figure or plt.figure(figsize=(10., 4.))
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
        p0 = X[i, :dyn.sv_th]
        d = np.array([np.cos(X[i, dyn.sv_th]), np.sin(X[i, dyn.sv_th])])
        p1, p2 = p0+P.l*d, p0-P.l*d
        _line_body.set_data([p1[0], p2[0]], [p1[1], p2[1]])
        time_text.set_text(time_template.format(i * dt))
        return time_text, _line_body
     
    dt = time[1]-time[0]
    dt_mili = dt*1000#25
    print(f'steps {len(time)} , {dt}, {dt_mili}')
    anim = animation.FuncAnimation(fig, animate, np.arange(1, len(time)),
                                   interval=dt_mili, blit=True, init_func=init, repeat_delay=200)

    return anim


def test_anim():
    fig = plt.figure()
    ax = plt.axes(xlim=(0, 4), ylim=(-2, 2))
    line, = ax.plot([], [], lw=3)

    def init():
        line.set_data([], [])
        return line,
    def animate(i):
        x = np.linspace(0, 4, 1000)
        y = np.sin(2 * np.pi * (x - 0.01 * i))
        line.set_data(x, y)
        return line,
    

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=200, interval=20, blit=True)
    return anim


def save_anim(filename, an, dt):
    print('encoding animation video, please wait, it will take a while')
    _start = time.time()
    fps = 1./dt/4; print(f'dt {dt} fps {fps}')
    an.save(filename, writer=animation.PillowWriter(fps=200)) # gif?
    _end = time.time()
    print(f'video encoded, saved to {filename}, Bye (took {_end-_start:.1f} s)')

    

import control
def test():
    P = dyn.Param()
    Xe, Ue = dyn.trim(P)
    A, B = dyn.jacobian(Xe, Ue, P)
    Q, R = [10, 1, 1, 1, 1, 1], [1, 1]
    (K, X, E) = control.lqr(A, B, np.diag(Q), np.diag(R))
    K = np.array(K) # control.lqr returns a np.matrix object, we use np.array
    cl_dyn = A - np.dot(B,K)
    print('cl poles: {}'.format(np.linalg.eig(cl_dyn)[0]))
    C = np.array([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0]])
    H = single_test1.get_precommand(A, B, C, K)
    X0 = [1, 0, 0, 0, 0, 0]
    time, X, U = single_test1.sim_state_feedback(P, X0, Ue, K, H, single_test1.step_x, tf=8.)
    Yc = [single_test1.step_x(_t) for _t in time]
    dyn.plot_trajectory(time, X, U, window_title="LQR")
    anim = animate(time, X, U, Yc, P)
    #anim = test_anim()
    #save_anim('/tmp/foo.apng', anim, 6/25.)
    plt.show()
  
if __name__ == "__main__":
    test()

