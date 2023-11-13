from BodySytem import BodySystem


def simulation_yarn():
    yarn_init_setup = {
        "name": "yarn",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.559064247131347, 0.559064247131347, -2 * 0.559064247131347],
        "vy": [0.349191558837891, 0.349191558837891, -2 * 0.349191558837891],
        "T": 55.5017624421301,
        "h": 0.0001,
    }

    yarn_euler = BodySystem(yarn_init_setup, ODESolver="Euler")
    yarn_euler.run_simulation()
    yarn_euler.plot_orbit()
    yarn_euler.plot_total_energy()
    yarn_euler.plot_angular_momentum()
    yarn_euler.plot_momentum_lineal()
    yarn_euler.save_error()

    yarn_rk4 = BodySystem(yarn_init_setup, ODESolver="RK4")
    yarn_rk4.run_simulation()
    yarn_rk4.plot_orbit()
    yarn_rk4.plot_total_energy()
    yarn_rk4.plot_angular_momentum()
    yarn_rk4.plot_momentum_lineal()
    yarn_rk4.save_error()

    yarn_rkf = BodySystem(yarn_init_setup, ODESolver="RKF")
    yarn_rkf.run_simulation()
    yarn_rkf.plot_orbit()
    yarn_rkf.plot_total_energy()
    yarn_rkf.plot_angular_momentum()
    yarn_rkf.plot_momentum_lineal()
    yarn_rkf.save_error()

    yarn_verlet = BodySystem(yarn_init_setup, ODESolver="Verlet")
    yarn_verlet.run_simulation()
    yarn_verlet.plot_orbit()
    yarn_verlet.plot_total_energy()
    yarn_verlet.plot_angular_momentum()
    yarn_verlet.plot_momentum_lineal()
    yarn_verlet.save_error()


if "__main__" == __name__:
    simulation_yarn()
