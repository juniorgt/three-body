from BodySytem import BodySystem


def simulation_goggles():
    goggles_init_setup = {
        "name": "goggles",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.0833000564575194, 0.0833000564575194, -2 * 0.0833000564575194],
        "vy": [0.127889282226563, 0.127889282226563, -2 * 0.127889282226563],
        "T": 10.4668176954385,
        "h":0.0001, # Mayor presicion para verlet
    }

    goggles_euler = BodySystem(goggles_init_setup, ODESolver="Euler")
    goggles_euler.run_simulation()
    goggles_euler.plot_orbit()
    goggles_euler.plot_total_energy()
    goggles_euler.plot_angular_momentum()
    goggles_euler.plot_momentum_lineal()
    goggles_euler.save_error()

    goggles_rk4 = BodySystem(goggles_init_setup, ODESolver="RK4")
    goggles_rk4.run_simulation()
    goggles_rk4.plot_orbit()
    goggles_rk4.plot_total_energy()
    goggles_rk4.plot_angular_momentum()
    goggles_rk4.plot_momentum_lineal()
    goggles_rk4.save_error()

    goggles_rkf = BodySystem(goggles_init_setup, ODESolver="RKF")
    goggles_rkf.run_simulation()
    goggles_rkf.plot_orbit()
    goggles_rkf.plot_total_energy()
    goggles_rkf.plot_angular_momentum()
    goggles_rkf.plot_momentum_lineal()
    goggles_rkf.save_error()

    goggles_verlet = BodySystem(goggles_init_setup, ODESolver="Verlet")
    goggles_verlet.run_simulation()
    goggles_verlet.plot_orbit()
    goggles_verlet.plot_total_energy()
    goggles_verlet.plot_angular_momentum()
    goggles_verlet.plot_momentum_lineal()
    goggles_verlet.save_error()


if "__main__" == __name__:
    simulation_goggles()
