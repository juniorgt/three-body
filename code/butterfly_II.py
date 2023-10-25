from BodySytem import BodySystem


def simulation_butterfly_II():
    butterflyII_init_setup = {
        "name": "Butterfly_II",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.392955223941802, 0.392955223941802, -2 * 0.392955223941802],
        "vy": [0.0975792352080344, 0.0975792352080344, -2 * 0.0975792352080344],
        "T": 7.00390738764014 * 2,
        "steps": 100000 * 2,
    }

    butterflyII_euler = BodySystem(butterflyII_init_setup, ODESolver="Euler")
    butterflyII_euler.run_simulation()
    butterflyII_euler.plot_orbit()
    butterflyII_euler.plot_total_energy()
    butterflyII_euler.plot_angular_momentum()
    butterflyII_euler.plot_momentum_lineal()
    butterflyII_euler.save_error()

    butterflyII_rk4 = BodySystem(butterflyII_init_setup, ODESolver="RK4")
    butterflyII_rk4.run_simulation()
    butterflyII_rk4.plot_orbit()
    butterflyII_rk4.plot_total_energy()
    butterflyII_rk4.plot_angular_momentum()
    butterflyII_rk4.plot_momentum_lineal()
    butterflyII_rk4.save_error()

    butterflyII_rkf = BodySystem(butterflyII_init_setup, ODESolver="RKF")
    butterflyII_rkf.run_simulation()
    butterflyII_rkf.plot_orbit()
    butterflyII_rkf.plot_total_energy()
    butterflyII_rkf.plot_angular_momentum()
    butterflyII_rkf.plot_momentum_lineal()
    butterflyII_rkf.save_error()

    butterflyII_verlet = BodySystem(butterflyII_init_setup, ODESolver="Verlet")
    butterflyII_verlet.run_simulation()
    butterflyII_verlet.plot_orbit()
    butterflyII_verlet.plot_total_energy()
    butterflyII_verlet.plot_angular_momentum()
    butterflyII_verlet.plot_momentum_lineal()
    butterflyII_verlet.save_error()


if "__main__" == __name__:
    simulation_butterfly_II()
