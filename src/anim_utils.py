import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import misc_utils as mu


def extrema(arr, d=0): return np.min(arr)-d, np.max(arr)+d
def subsample(l): return [_l[::4] if _l is not None else None for _l in l]  # FIXME hardcoded 100Hz -> 25fps 



class Animation:

    def __init__(self, time, X, U, P, Xsp=None, figure=None, ax=None):
        self.fig = figure or plt.figure(figsize=(10., 4.))
        time, X, U, Xsp =  subsample((time, X, U, Xsp))
        self.time, self.X, self.U, self.P, self.Xsp = time, X, U, P, Xsp # either that or pass it to animate?
        self.dt = time[1]-time[0]
        (_xmin, _xmax), (_ymin, _ymax) = self.compute_extends(X, P)
        self.ax = ax or self.fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(_xmin, _xmax),
                                             ylim=(_ymin, _ymax), facecolor=(0.5, 0.9, 0.9))

        self._time_template = 'time = {:0.1f}s'
        self._time_text = self.ax.text(0.025, 0.92, 'N/A', transform=self.ax.transAxes)


    def draw_quad(self, line, pos, theta, d):
        n = [d*np.cos(theta), d*np.sin(theta)]
        p1, p2 = pos + n, pos-n
        line.set_data([p1[0], p2[0]], [p1[1], p2[1]])
        
    def init(self):
        _lines = self.setup_drawing()
        return (self._time_text,) + _lines

    def animate(self, i):
        self._time_text.set_text(self._time_template.format(i * self.dt))
        _lines = self.animate_drawing(i)
        return (self._time_text,) + _lines
    
    def generate(self):
        dt_mili = self.dt*1000
        anim = animation.FuncAnimation(self.fig, self.animate, np.arange(1, len(self.time)),
                                       interval=dt_mili, blit=True, init_func=self.init, repeat_delay=200)
        return anim
