from utils.orbit_loader import get_orbit
from utils.simulations import run_simulations


def run_simulations_for_orbits(orbit_names: list, first_steps: list) -> None:
    for orbit_name in orbit_names:
        orbit = get_orbit(orbit_name)
        for h in first_steps:
            orbit.update({"h": h})
            run_simulations(orbit)


def main():
    orbit_names = [
        "Figure-8",
        "Butterfly_I",
        "Butterfly_II",
        "Bumblebee",
    ]

    first_steps = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]

    # run_simulations_for_orbits(orbit_names, first_steps)
    run_simulations_for_orbits(["Butterfly_II"], [1e-1, 1e-2])


if "__main__" == __name__:
    main()
