from BodySytem import BodySystem


def simulation_moth_II():
    moth_II_init_setup = {
        "name": "moth II",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.439165939331987, 0.439165939331987, -2 * 0.439165939331987],
        "vy": [0.452967645644678, 0.452967645644678, -2 * 0.452967645644678],
        "T": 28.6702783225658,
        "h": 0.001,
    }

    moth_II_euler = BodySystem(moth_II_init_setup, ODESolver="Euler")
    moth_II_euler.run_simulation()
    moth_II_euler.plot_orbit()
    moth_II_euler.plot_total_energy()
    moth_II_euler.plot_angular_momentum()
    moth_II_euler.plot_momentum_lineal()
    moth_II_euler.save_error()

    moth_II_rk4 = BodySystem(moth_II_init_setup, ODESolver="RK4")
    moth_II_rk4.run_simulation()
    moth_II_rk4.plot_orbit()
    moth_II_rk4.plot_total_energy()
    moth_II_rk4.plot_angular_momentum()
    moth_II_rk4.plot_momentum_lineal()
    moth_II_rk4.save_error()

    moth_II_rkf = BodySystem(moth_II_init_setup, ODESolver="RKF")
    moth_II_rkf.run_simulation()
    moth_II_rkf.plot_orbit()
    moth_II_rkf.plot_total_energy()
    moth_II_rkf.plot_angular_momentum()
    moth_II_rkf.plot_momentum_lineal()
    moth_II_rkf.save_error()

    moth_II_verlet = BodySystem(moth_II_init_setup, ODESolver="Verlet")
    moth_II_verlet.run_simulation()
    moth_II_verlet.plot_orbit()
    moth_II_verlet.plot_total_energy()
    moth_II_verlet.plot_angular_momentum()
    moth_II_verlet.plot_momentum_lineal()
    moth_II_verlet.save_error()


if "__main__" == __name__:
    simulation_moth_II()
