from BodySytem import BodySystem


def simulation_dragonfly():
    dragonfly_init_setup = {
        "name": "dragonfly",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.080584285736084, 0.080584285736084, -2 * 0.080584285736084],
        "vy": [0.588836087036132, 0.588836087036132, -2 * 0.588836087036132],
        "T": 21.2709751966648,
        "h":0.0001,
    }

    dragonfly_euler = BodySystem(dragonfly_init_setup, ODESolver="Euler")
    dragonfly_euler.run_simulation()
    dragonfly_euler.plot_orbit()
    dragonfly_euler.plot_total_energy()
    dragonfly_euler.plot_angular_momentum()
    dragonfly_euler.plot_momentum_lineal()
    dragonfly_euler.save_error()

    dragonfly_rk4 = BodySystem(dragonfly_init_setup, ODESolver="RK4")
    dragonfly_rk4.run_simulation()
    dragonfly_rk4.plot_orbit()
    dragonfly_rk4.plot_total_energy()
    dragonfly_rk4.plot_angular_momentum()
    dragonfly_rk4.plot_momentum_lineal()
    dragonfly_rk4.save_error()

    dragonfly_rkf = BodySystem(dragonfly_init_setup, ODESolver="RKF")
    dragonfly_rkf.run_simulation()
    dragonfly_rkf.plot_orbit()
    dragonfly_rkf.plot_total_energy()
    dragonfly_rkf.plot_angular_momentum()
    dragonfly_rkf.plot_momentum_lineal()
    dragonfly_rkf.save_error()

    dragonfly_verlet = BodySystem(dragonfly_init_setup, ODESolver="Verlet")
    dragonfly_verlet.run_simulation()
    dragonfly_verlet.plot_orbit()
    dragonfly_verlet.plot_total_energy()
    dragonfly_verlet.plot_angular_momentum()
    dragonfly_verlet.plot_momentum_lineal()
    dragonfly_verlet.save_error()


if "__main__" == __name__:
    simulation_dragonfly()
