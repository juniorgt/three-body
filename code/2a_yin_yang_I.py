from BodySytem import BodySystem


def simulation_yin_yang_I_2a():
    yin_yang_I_2a_init_setup = {
        "name": "yin_yang_I_2a",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.513938054919243, 0.513938054919243, -2 * 0.513938054919243],
        "vy": [0.304736003875733, 0.304736003875733, -2 * 0.304736003875733],
        "T": 17.328369755004,
        "h": 0.0001,
    }

    yin_yang_I_2a_euler = BodySystem(yin_yang_I_2a_init_setup, ODESolver="Euler")
    yin_yang_I_2a_euler.run_simulation()
    yin_yang_I_2a_euler.plot_orbit()
    yin_yang_I_2a_euler.plot_total_energy()
    yin_yang_I_2a_euler.plot_angular_momentum()
    yin_yang_I_2a_euler.plot_momentum_lineal()
    yin_yang_I_2a_euler.save_error()

    yin_yang_I_2a_rk4 = BodySystem(yin_yang_I_2a_init_setup, ODESolver="RK4")
    yin_yang_I_2a_rk4.run_simulation()
    yin_yang_I_2a_rk4.plot_orbit()
    yin_yang_I_2a_rk4.plot_total_energy()
    yin_yang_I_2a_rk4.plot_angular_momentum()
    yin_yang_I_2a_rk4.plot_momentum_lineal()
    yin_yang_I_2a_rk4.save_error()

    yin_yang_I_2a_rkf = BodySystem(yin_yang_I_2a_init_setup, ODESolver="RKF")
    yin_yang_I_2a_rkf.run_simulation()
    yin_yang_I_2a_rkf.plot_orbit()
    yin_yang_I_2a_rkf.plot_total_energy()
    yin_yang_I_2a_rkf.plot_angular_momentum()
    yin_yang_I_2a_rkf.plot_momentum_lineal()
    yin_yang_I_2a_rkf.save_error()

    yin_yang_I_2a_verlet = BodySystem(yin_yang_I_2a_init_setup, ODESolver="Verlet")
    yin_yang_I_2a_verlet.run_simulation()
    yin_yang_I_2a_verlet.plot_orbit()
    yin_yang_I_2a_verlet.plot_total_energy()
    yin_yang_I_2a_verlet.plot_angular_momentum()
    yin_yang_I_2a_verlet.plot_momentum_lineal()
    yin_yang_I_2a_verlet.save_error()


if "__main__" == __name__:
    simulation_yin_yang_I_2a()
