from utils.orbit_loader import get_orbit
from utils.simulations import run_simulation_rk4, run_simulations

orbit_names_test = [
    "Figure-8",
    "Butterfly-I",
]


first_steps_test = [1e-1, 1e-2]


def all_simulations():
    for orbit_name in orbit_names_test:
        orbit = get_orbit(orbit_name)
        for h in first_steps_test:
            orbit.update({"h": h})
            run_simulations(orbit)


run_simulation_rk4("Figure-8", h=1e-1)
