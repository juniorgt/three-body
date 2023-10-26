from BodySytem import BodySystem


def simulation_yin_yang_I_2b():
    yin_yang_I_2b_init_setup = {
        "name": "yin_yang_I_2b",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.282698682308198, 0.282698682308198, -2 * 0.282698682308198],
        "vy": [0.327208786129952, 0.327208786129952, -2 * 0.327208786129952],
        "T": 10.9625630756217,
        "h": 0.0001,
    }

    yin_yang_I_2b_euler = BodySystem(yin_yang_I_2b_init_setup, ODESolver="Euler")
    yin_yang_I_2b_euler.run_simulation()
    yin_yang_I_2b_euler.plot_orbit()
    yin_yang_I_2b_euler.plot_total_energy()
    yin_yang_I_2b_euler.plot_angular_momentum()
    yin_yang_I_2b_euler.plot_momentum_lineal()
    yin_yang_I_2b_euler.save_error()

    yin_yang_I_2b_rk4 = BodySystem(yin_yang_I_2b_init_setup, ODESolver="RK4")
    yin_yang_I_2b_rk4.run_simulation()
    yin_yang_I_2b_rk4.plot_orbit()
    yin_yang_I_2b_rk4.plot_total_energy()
    yin_yang_I_2b_rk4.plot_angular_momentum()
    yin_yang_I_2b_rk4.plot_momentum_lineal()
    yin_yang_I_2b_rk4.save_error()

    yin_yang_I_2b_rkf = BodySystem(yin_yang_I_2b_init_setup, ODESolver="RKF")
    yin_yang_I_2b_rkf.run_simulation()
    yin_yang_I_2b_rkf.plot_orbit()
    yin_yang_I_2b_rkf.plot_total_energy()
    yin_yang_I_2b_rkf.plot_angular_momentum()
    yin_yang_I_2b_rkf.plot_momentum_lineal()
    yin_yang_I_2b_rkf.save_error()

    yin_yang_I_2b_verlet = BodySystem(yin_yang_I_2b_init_setup, ODESolver="Verlet")
    yin_yang_I_2b_verlet.run_simulation()
    yin_yang_I_2b_verlet.plot_orbit()
    yin_yang_I_2b_verlet.plot_total_energy()
    yin_yang_I_2b_verlet.plot_angular_momentum()
    yin_yang_I_2b_verlet.plot_momentum_lineal()
    yin_yang_I_2b_verlet.save_error()


if "__main__" == __name__:
    simulation_yin_yang_I_2b()
