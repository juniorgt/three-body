from BodySytem import BodySystem


def simulation_butterfly_I():
    butterflyI_init_setup = {
        "name": "Butterfly_I",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.306892758965492, 0.306892758965492, -2 * 0.306892758965492],
        "vy": [0.125506782829762, 0.125506782829762, -2 * 0.125506782829762],
        "T": 6.23564136316479 * 1,
        "steps": 100000 * 1,
    }

    butterflyI_euler = BodySystem(butterflyI_init_setup, ODESolver="Euler")
    butterflyI_euler.run_simulation()
    butterflyI_euler.plot_orbit()
    butterflyI_euler.plot_total_energy()
    butterflyI_euler.plot_angular_momentum()
    butterflyI_euler.plot_momentum_lineal()
    butterflyI_euler.save_error()

    butterflyI_rk4 = BodySystem(butterflyI_init_setup, ODESolver="RK4")
    butterflyI_rk4.run_simulation()
    butterflyI_rk4.plot_orbit()
    butterflyI_rk4.plot_total_energy()
    butterflyI_rk4.plot_angular_momentum()
    butterflyI_rk4.plot_momentum_lineal()
    butterflyI_rk4.save_error()

    butterflyI_rkf = BodySystem(butterflyI_init_setup, ODESolver="RKF")
    butterflyI_rkf.run_simulation()
    butterflyI_rkf.plot_orbit()
    butterflyI_rkf.plot_total_energy()
    butterflyI_rkf.plot_angular_momentum()
    butterflyI_rkf.plot_momentum_lineal()
    butterflyI_rkf.save_error()

    butterflyI_verlet = BodySystem(butterflyI_init_setup, ODESolver="Verlet")
    butterflyI_verlet.run_simulation()
    butterflyI_verlet.plot_orbit()
    butterflyI_verlet.plot_total_energy()
    butterflyI_verlet.plot_angular_momentum()
    butterflyI_verlet.plot_momentum_lineal()
    butterflyI_verlet.save_error()


if "__main__" == __name__:
    simulation_butterfly_I()
