from bodySytem import BodySystem
from ODESolvers.methods import rk_methods, sym_methods
from utils.orbit_loader import get_orbit


def run_all_simulations(init_setup):
    for rk in rk_methods:
        BS = BodySystem(init_setup, ODESolver="rk", method=rk)
        BS.run_simulation()
        BS.plot()

    for sym in sym_methods:
        BS = BodySystem(init_setup, ODESolver="sym", method=sym)
        BS.run_simulation()
        BS.plot()


def run_simulation_rk4(orbit_name: str, h: float = None) -> None:
    orbit = get_orbit(orbit_name)
    if h is not None:
        orbit.update({"h": h})
    BS = BodySystem(orbit, ODESolver="rk", method="rk4")
    BS.run_simulation()
    BS.plot()


def run_simulation_fehlberg(orbit_name: str, h: float = None) -> None:
    orbit = get_orbit(orbit_name)
    if h is not None:
        orbit.update({"h": h})
    BS = BodySystem(orbit, ODESolver="rk", method="fehlberg")
    BS.run_simulation()
    BS.plot()


def run_simulation_dormand_prince(orbit_name: str, h: float = None) -> None:
    orbit = get_orbit(orbit_name)
    if h is not None:
        orbit.update({"h": h})
    BS = BodySystem(orbit, ODESolver="rk", method="dormand_prince")
    BS.run_simulation()
    BS.plot()


def run_simulation(
    orbit_name: str, ODESolver: str, method: str, h: float = None
) -> None:
    orbit = get_orbit(orbit_name)
    if h is not None:
        orbit.update({"h": h})
    BS = BodySystem(orbit, ODESolver, method)
    BS.run_simulation()
    BS.plot()
