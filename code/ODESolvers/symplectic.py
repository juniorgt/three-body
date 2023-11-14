import numpy as np


class SymIntegrator:
    def __init__(
        self,
        method,
        function,
        y0,
        t0=0.0,
        t_bound=300,
        first_step=1.0,
    ):
        self.method = method
        self.function = function
        self.y0 = y0
        self.t0 = t0
        self.first_step = first_step
        self.t_bound = t_bound

    def _single_step(self, t, y, h):
        for i in range(len(self.method[0])):
            y[::2] = y[::2] + h * np.array(y[1::2]) * self.method[0][i]
            y[1::2] = (
                y[1::2] + h * np.array(self.function(t, y)[1::2]) * self.method[1][i]
            )
        return t + h, y

    def run(self):
        stepsize = np.absolute(self.first_step)
        t = np.arange(self.t0, self.t_bound, self.first_step, dtype=float)
        y = np.empty((t.shape[0], len(self.y0)), dtype=float)
        y[0] = self.y0

        for i in range(1, t.shape[0]):
            t[i], y[i] = self._single_step(t[i - 1], y[i - 1], stepsize)
        return t, y
