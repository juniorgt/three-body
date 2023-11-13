from BodySytem import BodySystem


def simulation_moth_III():
    moth_III_init_setup = {
        "name": "moth_IIII",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.383443534851074, 0.383443534851074, -2 * 0.383443534851074],
        "vy": [0.377363693237305, 0.377363693237305, -2 * 0.377363693237305],
        "T": 25.8406180475758,
        "h": 0.001,
    }

    moth_III_euler = BodySystem(moth_III_init_setup, ODESolver="Euler")
    moth_III_euler.run_simulation()
    moth_III_euler.plot_orbit()
    moth_III_euler.plot_total_energy()
    moth_III_euler.plot_angular_momentum()
    moth_III_euler.plot_momentum_lineal()
    moth_III_euler.save_error()

    moth_III_rk4 = BodySystem(moth_III_init_setup, ODESolver="RK4")
    moth_III_rk4.run_simulation()
    moth_III_rk4.plot_orbit()
    moth_III_rk4.plot_total_energy()
    moth_III_rk4.plot_angular_momentum()
    moth_III_rk4.plot_momentum_lineal()
    moth_III_rk4.save_error()

    moth_III_rkf = BodySystem(moth_III_init_setup, ODESolver="RKF")
    moth_III_rkf.run_simulation()
    moth_III_rkf.plot_orbit()
    moth_III_rkf.plot_total_energy()
    moth_III_rkf.plot_angular_momentum()
    moth_III_rkf.plot_momentum_lineal()
    moth_III_rkf.save_error()

    moth_III_verlet = BodySystem(moth_III_init_setup, ODESolver="Verlet")
    moth_III_verlet.run_simulation()
    moth_III_verlet.plot_orbit()
    moth_III_verlet.plot_total_energy()
    moth_III_verlet.plot_angular_momentum()
    moth_III_verlet.plot_momentum_lineal()
    moth_III_verlet.save_error()


if "__main__" == __name__:
    simulation_moth_III()
