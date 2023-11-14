import numpy as np


class RKMethod:
    def __init__(
        self,
        method,
        function,
        y0,
        t0=0.0,
        t_bound=300,
        first_step=1.0,
        adaptive=False,
        tolerance=1.0e-20,
        endpoint=True,
    ):
        self.method = method
        self.function = function
        self.t0 = t0
        self.t_bound = t_bound
        self.y0 = y0
        self.first_step = first_step
        self.adaptive = adaptive
        self.tolerance = tolerance
        self.endpoint = endpoint

    def _step_impl(self, t, y, h, calc_error):
        """Perform a single Runge-Kutta step.

        Parameters
        ----------
        t : float
            Current time.
        y : ndarray, shape (n,)
            Current state.
        h : float
            Step to use.
        """

        k = np.zeros((len(y), self.method.shape[1] - 1))
        for i in range(self.method.shape[1] - 1):
            k[:, i] = self.function(
                t + h * self.method[i][0], y + h * k.dot(self.method[i, 1:])
            )

        y_new1 = y + h * k.dot(self.method[-1, 1:])
        if calc_error:
            y_new2 = y + h * k.dot(self.method[-2, 1:])
            error = np.absolute(y_new1 - y_new2)
            return t + h, y_new2, error
        else:
            return t + h, y_new1

    def run(self):
        y0 = np.asarray(self.y0)
        y0 = y0.astype(float, copy=False)
        if y0.ndim != 1:
            raise ValueError("'y0' must be 1-dimensional.")

        calc_error = True if self.method.shape[0] == self.method.shape[1] + 1 else False

        stepsize = np.absolute(self.first_step)

        t = [self.t0]
        y = [y0]
        if calc_error:
            error = [np.zeros(y0.shape[0])]

        while t[-1] < self.t_bound:
            if self.adaptive:
                stepsize = (
                    stepsize
                    * 0.9
                    * np.minimum(
                        np.maximum(
                            np.sqrt(
                                0.5
                                * np.amin(np.divide(self.tolerance, np.amax(error[-1])))
                            ),
                            0.3,
                        ),
                        2.0,
                    )
                )
            if self.endpoint:
                if t[-1] + stepsize > self.t_bound:
                    stepsize = np.absolute(self.t_bound - t[-1])

            h = stepsize
            new_values = self._step_impl(t[-1], y[-1], h, calc_error)
            t.append(new_values[0])
            y.append(new_values[1])
            if calc_error:
                error.append(new_values[2])

        y = np.array(y)
        return t, y


if __name__ == "__main__":

    def fun(t, y):
        return y - t**2 + 1

    y0 = 0.5

    rk = RKMethod(np.array([[0, 0], [None, 1.00000]]), fun, 0, 2, [y0], 0.2)
    a, b = rk.run()
    for x, y in zip(a, b):
        print(x, y)
