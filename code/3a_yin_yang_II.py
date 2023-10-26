from BodySytem import BodySystem


def simulation_yin_yang_II_3a():
    yin_yang_II_3a_init_setup = {
        "name": "yin_yang_II_3a",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.416822143554688, 0.416822143554688, -2 * 0.416822143554688],
        "vy": [0.330333312988282, 0.330333312988282, -2 * 0.330333312988282],
        "T":  55.78982856891 ,
        "h": 0.0001,
    }

    yin_yang_II_3a_euler = BodySystem(yin_yang_II_3a_init_setup, ODESolver="Euler")
    yin_yang_II_3a_euler.run_simulation()
    yin_yang_II_3a_euler.plot_orbit()
    yin_yang_II_3a_euler.plot_total_energy()
    yin_yang_II_3a_euler.plot_angular_momentum()
    yin_yang_II_3a_euler.plot_momentum_lineal()
    yin_yang_II_3a_euler.save_error()

    yin_yang_II_3a_rk4 = BodySystem(yin_yang_II_3a_init_setup, ODESolver="RK4")
    yin_yang_II_3a_rk4.run_simulation()
    yin_yang_II_3a_rk4.plot_orbit()
    yin_yang_II_3a_rk4.plot_total_energy()
    yin_yang_II_3a_rk4.plot_angular_momentum()
    yin_yang_II_3a_rk4.plot_momentum_lineal()
    yin_yang_II_3a_rk4.save_error()

    yin_yang_II_3a_rkf = BodySystem(yin_yang_II_3a_init_setup, ODESolver="RKF")
    yin_yang_II_3a_rkf.run_simulation()
    yin_yang_II_3a_rkf.plot_orbit()
    yin_yang_II_3a_rkf.plot_total_energy()
    yin_yang_II_3a_rkf.plot_angular_momentum()
    yin_yang_II_3a_rkf.plot_momentum_lineal()
    yin_yang_II_3a_rkf.save_error()

    yin_yang_II_3a_verlet = BodySystem(yin_yang_II_3a_init_setup, ODESolver="Verlet")
    yin_yang_II_3a_verlet.run_simulation()
    yin_yang_II_3a_verlet.plot_orbit()
    yin_yang_II_3a_verlet.plot_total_energy()
    yin_yang_II_3a_verlet.plot_angular_momentum()
    yin_yang_II_3a_verlet.plot_momentum_lineal()
    yin_yang_II_3a_verlet.save_error()


if "__main__" == __name__:
    simulation_yin_yang_II_3a()
