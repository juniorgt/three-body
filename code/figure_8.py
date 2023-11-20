from BodySytem import BodySystem


def simulation_figure_8():
    figure8_init_setup = {
        "name": "Figure-8",
        "ODESolver": "rk",
        "method": "original_rk",
        "G": 1,
        "M": [1, 1, 1],
        "y1": [-0.97000436, 0.4662036850, 0.24208753, 0.4323657300],  # x1, vx1, y1, vy1
        "y2": [0.0, -0.933240737, 0.0, -0.86473146],
        "y3": [0.97000436, 0.4662036850, -0.24208753, 0.4323657300],
        "T": 6.3259,
        "h": 5e-3,
    }

    figure8_euler = BodySystem(figure8_init_setup, ODESolver="rk", method="euler")
    figure8_euler.run_simulation()
    figure8_euler.plot_orbit()
    figure8_euler.plot_total_energy()
    figure8_euler.plot_angular_momentum()
    figure8_euler.plot_momentum_lineal()
    figure8_euler.save_error()

    figure8_rk4 = BodySystem(figure8_init_setup, ODESolver="rk", method="original_rk")
    figure8_rk4.run_simulation()
    figure8_rk4.plot_orbit()
    figure8_rk4.plot_total_energy()
    figure8_rk4.plot_angular_momentum()
    figure8_rk4.plot_momentum_lineal()
    figure8_rk4.save_error()

    figure8_rkf = BodySystem(figure8_init_setup, ODESolver="sym", method="euler")
    figure8_rkf.run_simulation()
    figure8_rkf.plot_orbit()
    figure8_rkf.plot_total_energy()
    figure8_rkf.plot_angular_momentum()
    figure8_rkf.plot_momentum_lineal()
    figure8_rkf.save_error()

    figure8_verlet = BodySystem(figure8_init_setup, ODESolver="sym", method="verlet")
    figure8_verlet.run_simulation()
    figure8_verlet.plot_orbit()
    figure8_verlet.plot_total_energy()
    figure8_verlet.plot_angular_momentum()
    figure8_verlet.plot_momentum_lineal()
    figure8_verlet.save_error()


if "__main__" == __name__:
    simulation_figure_8()
