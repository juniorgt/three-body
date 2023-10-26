from BodySytem import BodySystem


def simulation_yin_yang_II_3b():
    yin_yang_II_3b_init_setup = {
        "name": "yin_yang_II_3b",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.417342877101898, 0.417342877101898, -2 * 0.417342877101898],
        "vy": [0.313100116109848, 0.313100116109848, -2 * 0.313100116109848],
        "T": 54.2075992141846,
        "h": 0.0001,
    }

    yin_yang_II_3b_euler = BodySystem(yin_yang_II_3b_init_setup, ODESolver="Euler")
    yin_yang_II_3b_euler.run_simulation()
    yin_yang_II_3b_euler.plot_orbit()
    yin_yang_II_3b_euler.plot_total_energy()
    yin_yang_II_3b_euler.plot_angular_momentum()
    yin_yang_II_3b_euler.plot_momentum_lineal()
    yin_yang_II_3b_euler.save_error()

    yin_yang_II_3b_rk4 = BodySystem(yin_yang_II_3b_init_setup, ODESolver="RK4")
    yin_yang_II_3b_rk4.run_simulation()
    yin_yang_II_3b_rk4.plot_orbit()
    yin_yang_II_3b_rk4.plot_total_energy()
    yin_yang_II_3b_rk4.plot_angular_momentum()
    yin_yang_II_3b_rk4.plot_momentum_lineal()
    yin_yang_II_3b_rk4.save_error()

    yin_yang_II_3b_rkf = BodySystem(yin_yang_II_3b_init_setup, ODESolver="RKF")
    yin_yang_II_3b_rkf.run_simulation()
    yin_yang_II_3b_rkf.plot_orbit()
    yin_yang_II_3b_rkf.plot_total_energy()
    yin_yang_II_3b_rkf.plot_angular_momentum()
    yin_yang_II_3b_rkf.plot_momentum_lineal()
    yin_yang_II_3b_rkf.save_error()

    yin_yang_II_3b_verlet = BodySystem(yin_yang_II_3b_init_setup, ODESolver="Verlet")
    yin_yang_II_3b_verlet.run_simulation()
    yin_yang_II_3b_verlet.plot_orbit()
    yin_yang_II_3b_verlet.plot_total_energy()
    yin_yang_II_3b_verlet.plot_angular_momentum()
    yin_yang_II_3b_verlet.plot_momentum_lineal()
    yin_yang_II_3b_verlet.save_error()


if "__main__" == __name__:
    simulation_yin_yang_II_3b()
