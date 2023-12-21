from ODESolvers.methods import rk_methods
from ODESolvers.rungekutta import RKMethod


def fun(t, y):
    return y - t**2 + 1


y0 = 0.5
rk = RKMethod(rk_methods["euler"], fun, 0, 2, [y0], 0.2)
a, b = rk.run()
