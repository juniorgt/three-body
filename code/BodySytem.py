import os
import time
import uuid
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
import tomli_w
import tomllib
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

DEFAULT_SETUP = {
    "name": "Figure-8",
    "ODESolver": DEFAULT_ODE_SOLVER,
    "method": DEFAULT_METHOD,
    "G": 1,
    "M": [1.0, 1.0, 1.0],
    "y1": [-0.97000436, 0.4662036850, 0.24208753, 0.4323657300],  # x1, vx1, y1, vy1
    "y2": [0.0, -0.933240737, 0.0, -0.86473146],
    "y3": [0.97000436, 0.4662036850, -0.24208753, 0.4323657300],
    "T": 6.3259,
    "h": 5e-3,
}


class BodySystem:
    def __init__(
        self,
        init_setup=None,
        ODESolver: str | None = None,
        method: str | None = None,
    ):
        self._set_default_setup(init_setup)
        self._initialize_parameters(init_setup, ODESolver, method)
        self._create_save_directories()

    def _set_default_setup(self, init_setup):
        if init_setup is None:
            init_setup = DEFAULT_SETUP
        self.name = init_setup.get("name", f"orbit_{str(uuid.uuid4())[:5]}")
        self.G = init_setup.get("G", G)
        self.M = np.array(init_setup["M"], dtype=float)
        self.coords = np.concatenate(
            (init_setup["y1"], init_setup["y2"], init_setup["y3"]), dtype=float
        )
        self.T = init_setup.get("T", MAX_TIME)
        self.h = init_setup.get("h", DEFAULT_FIRST_STEP)
        self.steps = None

    def _initialize_parameters(
        self, init_setup, ODESolver: str | None, method: str | None
    ):
        self.ODESolver = (
            ODESolver
            if ODESolver is not None
            else init_setup.get("ODESolver", DEFAULT_ODE_SOLVER)
        )
        self.method = (
            method if method is not None else init_setup.get("method", DEFAULT_METHOD)
        )

    def _create_save_directories(self):
        self.save_route_images = os.path.join(
            save_route_images, self.name, f"{self.ODESolver}_{self.method}"
        )
        self.save_route_data = os.path.join(
            save_route_data, self.name, f"{self.ODESolver}_{self.method}"
        )
        os.makedirs(self.save_route_images, exist_ok=True)
        os.makedirs(self.save_route_data, exist_ok=True)

        self.path_t = os.path.join(self.save_route_data, "t.npy")
        self.path_y = os.path.join(self.save_route_data, "y.npy")

    def save_setup_to_toml(self):
        file_name = f"{self.name}_{self.ODESolver}_{self.method}.toml"
        self.file_path_toml = os.path.join(self.save_route_data, file_name)

        setup_data = {
            "init_setup": {
                "name": self.name,
                "ODESolver": self.ODESolver,
                "method": self.method,
                "G": self.G,
                "M": self.M.tolist(),
                "y1": self.coords[:4].tolist(),
                "y2": self.coords[4:8].tolist(),
                "y3": self.coords[8:].tolist(),
                "T": self.T,
                "h": self.h,
            }
        }

        try:
            with open(self.file_path_toml, "wb") as toml_file:
                tomli_w.dump(setup_data, toml_file)

        except Exception as e:
            print(f"Error saving the TOML file: {e}")

    def _generate_fun(
        self, masses: np.ndarray, G: float, nBodies=3
    ) -> Callable[[float, np.ndarray], np.ndarray]:
        n_dim = 2
        n_variable = 2

        def fun(t: float, y: np.ndarray) -> np.ndarray:
            body_states = y.reshape((nBodies, n_dim * n_variable)).astype(float)

            positions, velocities = body_states[:, ::2], body_states[:, 1::2]

            position_diff = positions[:, np.newaxis] - positions
            distances = np.linalg.norm(position_diff, axis=2)
            np.fill_diagonal(distances, 1)

            acc = (
                -G
                * masses[:, np.newaxis]
                * position_diff
                / (distances[:, :, np.newaxis] ** 3)
            )

            accelerations = np.sum(acc, axis=1)

            return np.vstack(
                (velocities.flatten(), accelerations.flatten())
            ).T.flatten()

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
        self.time, self.y = ODESolver.run()
        end_time = time.time()
        self.running_time = end_time - init_time
        self.save_setup_to_toml()
        self._save_running_time()
        self.save_simulation()
        self.maping_data(self.y)
        return self.time, self.y

    def maping_data(self, y):
        self.X = y[:, 0::4].T
        self.Y = y[:, 2::4].T
        self.VX = y[:, 1::4].T
        self.VY = y[:, 3::4].T

    def _save_running_time(self):
        try:
            with open(self.file_path_toml, "rb") as file:
                existing_data = tomllib.load(file)
        except FileNotFoundError:
            existing_data = {}

        metrics = {"metrics": {"running_time": self.running_time}}
        existing_data.update(metrics)

        try:
            with open(self.file_path_toml, "wb") as file:
                tomli_w.dump(existing_data, file)
        except Exception as e:
            print(f"Error saving the TOML file: {e}")

    def save_simulation(self) -> None:
        try:
            with open(self.path_t, "wb") as file_time:
                np.save(file_time, self.time)
        except Exception as e:
            print(f"Error when saving time: {e}")

        try:
            with open(self.path_y, "wb") as file_y:
                np.save(file_y, self.y)
        except Exception as e:
            print(f"Error saving 'y' data: {e}")

    def load_simulation(self):
        if not (os.path.exists(self.path_t) and os.path.exists(self.path_y)):
            raise FileNotFoundError("The specified files do not exist!")

        try:
            with open(self.path_t, "rb") as f:
                self.time = np.load(f)

            with open(self.path_y, "rb") as f:
                self.y = np.load(f, allow_pickle=True)

        except Exception as e:
            print(f"Error loading files: {e}")

        self.maping_data(self.y)

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
        plt.close(fig)

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
        with open(self.file_path_toml, "rb") as f:
            data = tomllib.load(f)
        dic_metrics = {
            "metrics": {
                "running_time": self.running_time,
                "error_enery": self._error(self.total_energy),
                "error_angular_momentum": self._error(self.angular_momentum),
            }
        }
        data |= dic_metrics
        with open(self.file_path_toml, "wb") as f:
            tomli_w.dump(data, f)

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
            os.path.join(
                route,
                f"{self.name}_{self.ODESolver}_{self.method}_angular_momentum.png",
            )
        )
        plt.close(fig)

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
            os.path.join(
                route, f"{self.name}_{self.ODESolver}_{self.method}_total_energy.png"
            )
        )
        plt.close(fig)

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
            os.path.join(
                route, f"{self.name}_{self.ODESolver}_{self.method}_linear_momentum.png"
            )
        )
        plt.close(fig)

    def plot(self):
        self.plot_orbit()
        self.plot_total_energy()
        self.plot_angular_momentum()
        self.plot_momentum_lineal()
        self.save_error()


if "__main__" == __name__:
    init_setup = {
        "name": "Figure-8",
        "ODESolver": DEFAULT_ODE_SOLVER,
        "G": 1,
        "M": [1, 1, 1],
        "y1": [-0.97000436, 0.4662036850, 0.24208753, 0.4323657300],  # x1, vx1, y1, vy1
        "y2": [0.0, -0.933240737, 0.0, -0.86473146],
        "y3": [0.97000436, 0.4662036850, -0.24208753, 0.4323657300],
        "T": 6.3259,
        "h": 5e-5,
    }

    bdrk = BodySystem(init_setup=init_setup, ODESolver="rk", method="four")
    # t, y = bdrk.run_simulation()
    # bdrk.save_simulation()
    # bdrk.load_simulation()
    # print(bdrk.time)
    # print(bdrk.y)
    # bdrk.plot_orbit()
    # bdsym = BodySystem(ODESolver="sym", method="verlet")
    # t, y = bdsym.run_simulation()
    # bdsym.plot_orbit()
