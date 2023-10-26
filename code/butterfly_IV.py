from BodySytem import BodySystem


def simulation_butterfly_IV():
    butterfly_IV_init_setup = {
        "name": "butterfly_IV",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-1.0, 1.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "vx": [0.350112121391296, 0.350112121391296, -2 * 0.350112121391296],
        "vy": [0.0793394773483276, 0.0793394773483276, -2 * 0.0793394773483276],
        "T": 79.4758748952101 ,
        "h":0.0001, # Mucha mayor presiccion > 0.0001
    }

    butterfly_IV_euler = BodySystem(butterfly_IV_init_setup, ODESolver="Euler")
    butterfly_IV_euler.run_simulation()
    butterfly_IV_euler.plot_orbit()
    butterfly_IV_euler.plot_total_energy()
    butterfly_IV_euler.plot_angular_momentum()
    butterfly_IV_euler.plot_momentum_lineal()
    butterfly_IV_euler.save_error()

    butterfly_IV_rk4 = BodySystem(butterfly_IV_init_setup, ODESolver="RK4")
    butterfly_IV_rk4.run_simulation()
    butterfly_IV_rk4.plot_orbit()
    butterfly_IV_rk4.plot_total_energy()
    butterfly_IV_rk4.plot_angular_momentum()
    butterfly_IV_rk4.plot_momentum_lineal()
    butterfly_IV_rk4.save_error()

    butterfly_IV_rkf = BodySystem(butterfly_IV_init_setup, ODESolver="RKF")
    butterfly_IV_rkf.run_simulation()
    butterfly_IV_rkf.plot_orbit()
    butterfly_IV_rkf.plot_total_energy()
    butterfly_IV_rkf.plot_angular_momentum()
    butterfly_IV_rkf.plot_momentum_lineal()
    butterfly_IV_rkf.save_error()

    butterfly_IV_verlet = BodySystem(butterfly_IV_init_setup, ODESolver="Verlet")
    butterfly_IV_verlet.run_simulation()
    butterfly_IV_verlet.plot_orbit()
    butterfly_IV_verlet.plot_total_energy()
    butterfly_IV_verlet.plot_angular_momentum()
    butterfly_IV_verlet.plot_momentum_lineal()
    butterfly_IV_verlet.save_error()


if "__main__" == __name__:
    simulation_butterfly_IV()
