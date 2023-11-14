import json
import os
import time
import uuid

import matplotlib.pyplot as plt
import numpy as np
from ODESolvers.methods import rk_methods, sym_methods
from ODESolvers.rungekutta import RKMethod
from ODESolvers.symplectic import SymIntegrator
from routes import save_route_data, save_route_images

EPS = np.finfo(float).eps

G = 6.67408313131313e-11
MAX_TIME = 300
DEFAULT_ODE_SOLVER = "rk"
DEFAULT_METHOD = "original_rk"
DEFAULT_FIRST_STEP = 5e-3

default_setup = {
    "name": "Figure-8",
    "ODESolver": DEFAULT_ODE_SOLVER,
    "G": 1,
    "M": [1, 1, 1],
    "x": [-0.97000436, 0.0, 0.97000436],
    "y": [0.24208753, 0.0, -0.24208753],
    "vx": [0.4662036850, -0.933240737, 0.4662036850],
    "vy": [0.4323657300, -0.86473146, 0.4323657300],
    "T": 6.3259,
    "h": 5e-3,
}


class BodySystem:
    def __init__(self, init_setup=None, ODESolver=None, method=None):
        if init_setup is None:
            init_setup = default_setup

        self.name = init_setup.get("name", f"orbit_{str(uuid.uuid4())[:5]}")

        self.ODESolver = (
            ODESolver
            if ODESolver is not None
            else init_setup.get("ODESolver", DEFAULT_ODE_SOLVER)
        )

        self.method = (
            method if method is not None else init_setup.get("method", DEFAULT_METHOD)
        )

        self.G = init_setup.get("G", 6.67408313131313e-11)
        self.M = init_setup["M"]
        self.x = init_setup["x"]
        self.y = init_setup["y"]
        self.vx = init_setup["vx"]
        self.vy = init_setup["vy"]
        self.coords = np.concatenate((self.x, self.y, self.vx, self.vy), dtype=float)

        self.T = init_setup.get("T", MAX_TIME)
        self.h = init_setup.get("h", DEFAULT_FIRST_STEP)
        self.steps = None

        self.save_route_images = os.path.join(save_route_images, self.name)
        self.save_route_data = os.path.join(save_route_data, self.name)
        self._make_dir_saves()
        self._encode_init_setup()
        self.running_time = 0

    def _make_dir_saves(self):
        os.makedirs(self.save_route_images, exist_ok=True)
        os.makedirs(self.save_route_data, exist_ok=True)

    def _encode_init_setup(self):
        init_setup = {
            "ODESolver": self.ODESolver,
            "G": self.G,
            "m1": self.M[0],
            "v1": [self.vx[0], self.vy[0]],
            "r1": [self.x[0], self.y[0]],
            "m2": self.M[1],
            "v2": [self.vx[1], self.vy[1]],
            "r2": [self.x[1], self.y[1]],
            "m3": self.M[2],
            "v3": [self.vx[2], self.vy[2]],
            "r3": [self.x[2], self.y[2]],
            "T": self.T,
            "steps": self.steps,
            "h": self.h,
        }
        json_init_setup = json.dumps(init_setup, indent=4)
        # with open(
        #     os.path.join(self.save_route_images, f"{self.name}_{self.ODESolver}.json"),
        #     "w",
        # ) as json_file:
        #     json_file.write(json_init_setup)

        with open(
            os.path.join(self.save_route_data, f"{self.name}_{self.ODESolver}.json"),
            "w",
        ) as json_file:
            json_file.write(json_init_setup)

    "TODO: Correguir para que funcione con los metodo symplecticos"

    def _generate_fun(self, masses, G, nBodies=3):
        def fun(t, y):
            rx, ry, vx, vy = (
                y[:nBodies],
                y[nBodies : 2 * nBodies],
                y[2 * nBodies : 3 * nBodies],
                y[3 * nBodies :],
            )
            acc = np.zeros_like(y)
            for n in range(nBodies):
                xn, yn = rx[n], ry[n]
                acc_vx, acc_vy = 0.0, 0.0
                for i in range(nBodies):
                    if i != n:
                        sep = np.sqrt((xn - rx[i]) ** 2 + (yn - ry[i]) ** 2)
                        acc_vx -= G * masses[i] * (xn - rx[i]) / sep**3
                        acc_vy -= G * masses[i] * (yn - ry[i]) / sep**3
                acc[n] = vx[n]
                acc[n + nBodies] = vy[n]
                acc[n + 2 * nBodies] = acc_vx
                acc[n + 3 * nBodies] = acc_vy

            return acc

        return fun

    def run_simulation(self):
        fun = self._generate_fun(self.M, self.G, nBodies=3)
        if self.ODESolver == "rk":
            ODESolver = RKMethod(
                rk_methods[self.method],
                fun,
                self.coords,
                0,
                self.T,
                self.h,
            )
        elif self.ODESolver == "sym":
            ODESolver = SymIntegrator(
                sym_methods[self.method], fun, self.coords, 0, self.T, self.h
            )
        else:
            raise ValueError("No se especifico ODESolver")
        init_time = time.time()
        self.time, y = ODESolver.run()
        end_time = time.time()
        # print(self.time)
        # print(y)
        self.running_time = end_time - init_time
        self._save_runnig_time()
        self.running_time = end_time - init_time
        self.X = y[:, :3].T
        self.Y = y[:, 3:6].T
        self.VX = y[:, 6:9].T
        self.VY = y[:, 9:12].T
        return self.time, self.X, self.Y, self.VX, self.VY

    def _save_runnig_time(self):
        save_path = os.path.join(
            self.save_route_data, f"{self.name}_{self.ODESolver}_runing_time.txt"
        )
        with open(save_path, "w") as file:
            file.write(f"Running time = {self.running_time}\n")

    def save_simulation(self):
        pass

    def load_simulation(self):
        pass

    def _create_plot_orbit(self):
        fig, ax = plt.subplots()
        fig.suptitle(f"{self.name}")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.plot(self.X[0, :], self.Y[0, :], color="b", label="cuerpo 1")
        ax.plot(self.X[1, :], self.Y[1, :], color="r", label="cuerpo 2")
        ax.plot(self.X[2, :], self.Y[2, :], color="k", label="cuerpo 3")
        ax.scatter(
            self.X[:, -1],
            self.Y[:, -1],
            marker="o",
            s=100,
            edgecolor="none",
            c=["b", "r", "k"],
        )
        # Add a legend
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height])
        ax.legend(loc="center right", bbox_to_anchor=(1.35, 0.5))
        return fig, ax

    def plot_orbit(self, route=None):
        if route is None:
            route = self.save_route_images
        fig, _ = self._create_plot_orbit()
        fig.savefig(
            os.path.join(route, f"{self.name}_{self.ODESolver}_{self.method}.png")
        )

    def cal_angular_momentum(self):
        angular_momentum = np.zeros(len(self.time))
        for i in range(len(self.time)):
            for j in range(3):
                mi = self.M[j]
                xi, yi = self.X[j, i], self.Y[j, i]
                vxi, vyi = self.VX[j, i], self.VY[j, i]
                angular_momentum[i] += mi * (xi * vyi - yi * vxi)
        return angular_momentum

    def _error(self, arr):
        i = 0
        a = arr[i]
        while a == 0:
            a = arr[i + 1]
            i += 1
        b = np.abs(arr[1:] - a / a)
        return np.mean(b)

    def save_error(self):
        save_path = os.path.join(
            self.save_route_data, f"{self.name}_{self.ODESolver}_runing_time.txt"
        )
        with open(save_path, "w") as file:
            file.write(f"Running time = {self.running_time}\n")
            file.write(f"Error Enery = {self._error(self.total_energy)}\n")
            file.write(
                f"Error momento angular = {self._error(self.angular_momentum)}\n"
            )

    def _create_plot_angular_momentum(self):
        self.angular_momentum = self.cal_angular_momentum()
        fig, ax = plt.subplots()
        fig.suptitle("Tiempo vs Momento Angular")
        ax.plot(self.time, self.angular_momentum)
        ax.set_xlabel("Tiempo")
        ax.set_ylabel("Momento angular")
        return fig, ax

    def plot_angular_momentum(self, route=None):
        if route is None:
            route = self.save_route_images
        fig, _ = self._create_plot_angular_momentum()
        fig.savefig(
            os.path.join(route, f"{self.name}_{self.ODESolver}_angular_momentum.png")
        )

    def cal_total_energy(self):
        kinetic_energy = np.zeros(len(self.time))
        potential_energy = np.zeros(len(self.time))
        for i in range(len(self.time)):
            for j in range(3):
                mi = self.M[j]
                vxi, vyi = self.VX[j, i], self.VY[j, i]
                xi, yi = self.X[j, i], self.Y[j, i]
                kinetic_energy[i] += 0.5 * mi * (vxi**2 + vyi**2)
                for k in range(3):
                    if k != j:
                        mj = self.M[k]
                        xj, yj = self.X[k, i], self.Y[k, i]
                        r = np.sqrt((xi - xj) ** 2 + (yi - yj) ** 2)
                        potential_energy[i] += -G * mi * mj / r

        total_energy = kinetic_energy + potential_energy
        return total_energy

    def _create_plot_total_energy(self):
        self.total_energy = self.cal_total_energy()
        fig, ax = plt.subplots()
        fig.suptitle("Tiempo vs Energia Total")
        ax.plot(self.time, self.total_energy)
        ax.set_xlabel("Tiempo")
        ax.set_ylabel("Energia Total")
        return fig, ax

    def plot_total_energy(self, route=None):
        if route is None:
            route = self.save_route_images
        fig, _ = self._create_plot_total_energy()
        fig.savefig(
            os.path.join(route, f"{self.name}_{self.ODESolver}_total_energy.png")
        )

    def cal_linear_momentum(self):
        linear_momentum_x = np.zeros(len(self.time))
        linear_momentum_y = np.zeros(len(self.time))

        for i in range(len(self.time)):
            for j in range(3):
                mi = self.M[j]
                vxi, vyi = self.VX[j, i], self.VY[j, i]
                linear_momentum_x[i] += mi * vxi
                linear_momentum_y[i] += mi * vyi

        return linear_momentum_x, linear_momentum_y

    def _create_plot_linear_momentum(self):
        linear_momentum_x, linear_momentum_y = self.cal_linear_momentum()
        fig, ax = plt.subplots(2, 1)
        fig.suptitle("Tiempo vs Momento Lineal")
        ax[0].plot(self.time, linear_momentum_x, label="Linear Momentum X", color="b")
        ax[1].plot(self.time, linear_momentum_y, label="Linear Momentum Y", color="r")
        ax[0].set_xlabel("Tiempo")
        ax[1].set_xlabel("Tiempo")
        ax[0].set_ylabel("Momento Lineal")
        ax[1].set_ylabel("Momento Lineal")
        ax[0].legend()
        ax[1].legend()
        plt.tight_layout()
        return fig, ax

    def _create_plot_linear_momentumv2(self):
        linear_momentum_x, linear_momentum_y = self.cal_linear_momentum()
        fig, ax = plt.subplots()
        ax.plot(self.time, linear_momentum_x, label="Momento Lineal X")
        ax.plot(self.time, linear_momentum_y, label="Momento Lineal Y")
        ax.set_yscale("log")
        ax.set_xlabel("Tiempo")
        ax.set_ylabel("Momento Lineal")
        ax.legend()
        return fig, ax

    def plot_momentum_lineal(self, route=None):
        if route is None:
            route = self.save_route_images
        fig, _ = self._create_plot_linear_momentum()
        fig.savefig(
            os.path.join(route, f"{self.name}_{self.ODESolver}_linear_momentum.png")
        )


if "__main__" == __name__:
    init_setup = {
        "name": "Figure-8",
        "G": 1,
        "M": [1, 1, 1],
        "x": [-0.97000436, 0.0, 0.97000436],
        "y": [0.24208753, 0.0, -0.24208753],
        "vx": [0.4662036850, -0.933240737, 0.4662036850],
        "vy": [0.4323657300, -0.86473146, 0.4323657300],
        "T": 6.3259,
        "h": 2,
    }

    bdrk = BodySystem(init_setup, ODESolver="rk", method="euler")
    t, x, y, vx, vy = bdrk.run_simulation()
    print(t)
    print(x)
    bdrk.plot_orbit()
    bdsym = BodySystem(init_setup, ODESolver="sym", method="euler")
    t, x, y, vx, vy = bdsym.run_simulation()
    print(t)
    print(x)
    bdrk.plot_orbit()
    # bs1.plot_angular_momentum()
    # bs1.plot_total_energy()
    # bs1.plot_momentum_lineal()
