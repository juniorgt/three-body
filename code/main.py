from utils.orbit_loader import get_orbit
from utils.simulations import run_all_simulations, run_simulation


def run_simulations_for_orbits(orbit_names: list, first_steps: list) -> None:
    for orbit_name in orbit_names:
        orbit = get_orbit(orbit_name)
        for h in first_steps:
            orbit.update({"h": h})
            run_all_simulations(orbit)


def main():
    orbit_names = [
        "Figure-8",
        "Butterfly-I",
        "Butterfly-II",
        "Bumblebee",
        "Moth-I",
        "Moth-II",
        "Butterfly-III",
        "Goggles",
        # "Butterfly-IV",
        "Dragonfly",
        "Yarn",
        "2a-Yin-yang-I",
        "2b-Yin-yang-I",
        # "3a-Yin-yang-II",
        # "3b-Yin-yang-II",
    ]

    first_steps = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]

    run_simulations_for_orbits(orbit_names, first_steps)
    # run_simulations_for_orbits(["3b-Yin-yang-II"], [1e-5])


if "__main__" == __name__:
    # main()
    run_simulation("3b-Yin-yang-II", "rk", "dormand_prince", 1e-4)
