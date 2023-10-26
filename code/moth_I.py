from BodySytem import BodySystem


def simulation_moth_I():
    moth_I_init_setup = {
        "name": "moth_I",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.464445237398184, 0.464445237398184, -2 * 0.464445237398184],
        "vy": [0.396059973403921, 0.396059973403921, -2 * 0.396059973403921],
        "T": 14.8939113169584,
        "h": 0.001,
    }

    moth_I_euler = BodySystem(moth_I_init_setup, ODESolver="Euler")
    moth_I_euler.run_simulation()
    moth_I_euler.plot_orbit()
    moth_I_euler.plot_total_energy()
    moth_I_euler.plot_angular_momentum()
    moth_I_euler.plot_momentum_lineal()
    moth_I_euler.save_error()

    moth_I_rk4 = BodySystem(moth_I_init_setup, ODESolver="RK4")
    moth_I_rk4.run_simulation()
    moth_I_rk4.plot_orbit()
    moth_I_rk4.plot_total_energy()
    moth_I_rk4.plot_angular_momentum()
    moth_I_rk4.plot_momentum_lineal()
    moth_I_rk4.save_error()

    moth_I_rkf = BodySystem(moth_I_init_setup, ODESolver="RKF")
    moth_I_rkf.run_simulation()
    moth_I_rkf.plot_orbit()
    moth_I_rkf.plot_total_energy()
    moth_I_rkf.plot_angular_momentum()
    moth_I_rkf.plot_momentum_lineal()
    moth_I_rkf.save_error()

    moth_I_verlet = BodySystem(moth_I_init_setup, ODESolver="Verlet")
    moth_I_verlet.run_simulation()
    moth_I_verlet.plot_orbit()
    moth_I_verlet.plot_total_energy()
    moth_I_verlet.plot_angular_momentum()
    moth_I_verlet.plot_momentum_lineal()
    moth_I_verlet.save_error()


if "__main__" == __name__:
    simulation_moth_I()
