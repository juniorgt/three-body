from BodySytem import BodySystem


def simulation_bumblebee():
    bumblebee_init_setup = {
        "name": "bumblebee",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.1842784887, 0.1842784887, -2 * 0.1842784887],
        "vy": [0.5871881740, 0.5871881740, -2 * 0.5871881740],
        "T": 63.5343529785,
        "steps": 1000000,
    }

    bumblebee_euler = BodySystem(bumblebee_init_setup, ODESolver="Euler")
    bumblebee_euler.run_simulation()
    bumblebee_euler.plot_orbit()
    bumblebee_euler.plot_total_energy()
    bumblebee_euler.plot_angular_momentum()
    bumblebee_euler.plot_momentum_lineal()
    bumblebee_euler.save_error()

    bumblebee_rk4 = BodySystem(bumblebee_init_setup, ODESolver="RK4")
    bumblebee_rk4.run_simulation()
    bumblebee_rk4.plot_orbit()
    bumblebee_rk4.plot_total_energy()
    bumblebee_rk4.plot_angular_momentum()
    bumblebee_rk4.plot_momentum_lineal()
    bumblebee_rk4.save_error()

    bumblebee_rkf = BodySystem(bumblebee_init_setup, ODESolver="RKF")
    bumblebee_rkf.run_simulation()
    bumblebee_rkf.plot_orbit()
    bumblebee_rkf.plot_total_energy()
    bumblebee_rkf.plot_angular_momentum()
    bumblebee_rkf.plot_momentum_lineal()
    bumblebee_rkf.save_error()

    bumblebee_verlet = BodySystem(bumblebee_init_setup, ODESolver="Verlet")
    bumblebee_verlet.run_simulation()
    bumblebee_verlet.plot_orbit()
    bumblebee_verlet.plot_total_energy()
    bumblebee_verlet.plot_angular_momentum()
    bumblebee_verlet.plot_momentum_lineal()
    bumblebee_verlet.save_error()


if "__main__" == __name__:
    simulation_bumblebee()
