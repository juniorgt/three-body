from BodySytem import BodySystem


def simulation_figure_8():
    figure8_init_setup = {
        "name": "Figure-8_2",
        "ODESolver": "Euler",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-0.97000436, 0.0, 0.97000436],
        "y": [0.24208753, 0.0, -0.24208753],
        "vx": [0.4662036850, -0.933240737, 0.4662036850],
        "vy": [0.4323657300, -0.86473146, 0.4323657300],
        "T": 6.3259 * 200,
        "steps": 1000 * 200,
    }

    figure8_euler = BodySystem(figure8_init_setup, ODESolver="Euler")
    figure8_euler.run_simulation()
    figure8_euler.plot_orbit()
    figure8_euler.plot_total_energy()
    figure8_euler.plot_angular_momentum()
    figure8_euler.plot_momentum_lineal()

    figure8_rk4 = BodySystem(figure8_init_setup, ODESolver="RK4")
    figure8_rk4.run_simulation()
    figure8_rk4.plot_orbit()
    figure8_rk4.plot_total_energy()
    figure8_rk4.plot_angular_momentum()
    figure8_rk4.plot_momentum_lineal()

    figure8_rk5 = BodySystem(figure8_init_setup, ODESolver="RK5")
    figure8_rk5.run_simulation()
    figure8_rk5.plot_orbit()
    figure8_rk5.plot_total_energy()
    figure8_rk5.plot_angular_momentum()
    figure8_rk5.plot_momentum_lineal()

    figure8_verlet = BodySystem(figure8_init_setup, ODESolver="Verlet")
    figure8_verlet.run_simulation()
    figure8_verlet.plot_orbit()
    figure8_verlet.plot_total_energy()
    figure8_verlet.plot_angular_momentum()
    figure8_verlet.plot_momentum_lineal()


if "__main__" == __name__:
    simulation_figure_8()
