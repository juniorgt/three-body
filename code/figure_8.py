from bodySytem import BodySystem
from ODESolvers.methods import rk_methods, sym_methods


def simulation_figure_8(h):
    init_setup = {
        "name": "Figure-8",
        "ODESolver": "rk",
        "method": "original_rk",
        "G": 1,
        "M": [1, 1, 1],
        "y1": [-0.97000436, 0.4662036850, 0.24208753, 0.4323657300],
        "y2": [0.0, -0.933240737, 0.0, -0.86473146],
        "y3": [0.97000436, 0.4662036850, -0.24208753, 0.4323657300],
        "T": 6.3259,
        "h": h,
    }

    for rk in rk_methods.keys():
        BS = BodySystem(init_setup, ODESolver="rk", method=rk)
        BS.run_simulation()
        BS.plot()

    for sym in sym_methods.keys():
        BS = BodySystem(init_setup, ODESolver="sym", method=sym)
        BS.run_simulation()
        BS.plot()


if "__main__" == __name__:
    list_first_steps = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
    for h in list_first_steps:
        simulation_figure_8(h)
