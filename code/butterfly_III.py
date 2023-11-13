from BodySytem import BodySystem


def simulation_butterfly_III():
    butterfly_III_init_setup = {
        "name": "butterfly_III",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.405915588857606, 0.405915588857606, -2 * 0.405915588857606],
        "vy": [0.452967645644678, 0.452967645644678, -2 * 0.452967645644678],
        "T": 13.8657626785699,
        "h": 0.0001,  # Mayor presicion
    }

    butterfly_III_euler = BodySystem(butterfly_III_init_setup, ODESolver="Euler")
    butterfly_III_euler.run_simulation()
    butterfly_III_euler.plot_orbit()
    butterfly_III_euler.plot_total_energy()
    butterfly_III_euler.plot_angular_momentum()
    butterfly_III_euler.plot_momentum_lineal()
    butterfly_III_euler.save_error()

    butterfly_III_rk4 = BodySystem(butterfly_III_init_setup, ODESolver="RK4")
    butterfly_III_rk4.run_simulation()
    butterfly_III_rk4.plot_orbit()
    butterfly_III_rk4.plot_total_energy()
    butterfly_III_rk4.plot_angular_momentum()
    butterfly_III_rk4.plot_momentum_lineal()
    butterfly_III_rk4.save_error()

    butterfly_III_rkf = BodySystem(butterfly_III_init_setup, ODESolver="RKF")
    butterfly_III_rkf.run_simulation()
    butterfly_III_rkf.plot_orbit()
    butterfly_III_rkf.plot_total_energy()
    butterfly_III_rkf.plot_angular_momentum()
    butterfly_III_rkf.plot_momentum_lineal()
    butterfly_III_rkf.save_error()

    butterfly_III_verlet = BodySystem(butterfly_III_init_setup, ODESolver="Verlet")
    butterfly_III_verlet.run_simulation()
    butterfly_III_verlet.plot_orbit()
    butterfly_III_verlet.plot_total_energy()
    butterfly_III_verlet.plot_angular_momentum()
    butterfly_III_verlet.plot_momentum_lineal()
    butterfly_III_verlet.save_error()


if "__main__" == __name__:
    simulation_butterfly_III()
