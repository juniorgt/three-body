from bodySytem import BodySystem
from ODESolvers.methods import rk_methods, sym_methods


def simulation_butterfly_I(h):
    init_setup = {
        "name": "Butterfly-I",
        "G": 1,
        "M": [1, 1, 1],
        "y1": [-1.0, 0.306892758965492, 0.0, 0.125506782829762],
        "y2": [1.0, 0.306892758965492, 0.0, 0.125506782829762],
        "y3": [0.0, -2 * 0.306892758965492, 0.0, -2 * 0.125506782829762],
        "T": 6.23564136316479,
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
        simulation_butterfly_I(h)
