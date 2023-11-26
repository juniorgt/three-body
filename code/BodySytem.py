import os
import time
import uuid
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
import tomli_w
import tomllib
from matplotlib.colors import to_hex
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

        file_name = f"{self.name}_{self.ODESolver}_{self.method}.text"
        self.path_data = os.path.join(self.save_route_data, file_name)

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
        self._maping_data(self.y)
        self.calculate_total_energy()
        return self.time, self.y

    def _maping_data(self, y):
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
        simulation_data = np.hstack((self.time[:, np.newaxis], self.y))
        try:
            np.savetxt(self.path_data, simulation_data)
        except OSError:
            print(f"Error saving data to {self.path_data}")

    def load_simulation(self) -> None:
        if not os.path.exists(self.path_data):
            raise FileNotFoundError("The specified file does not exist!")

        try:
            self.data = np.loadtxt(self.path_data)
            self.time = self.data[:, 0]
            self.y = self.data[:, 1:]
            self._maping_data(self.y)
        except Exception as e:
            print(f"Error loading file: {e}")

    def _create_plot_orbit(
        self,
        marker="o",
        s=200,
        figsize=(9.6, 7.2),
        dpi=200,
        fontsize=14,
        legend_loc="center right",
        legend_bbox_to_anchor=(1.35, 0.5),
        linewidth=3,
    ):
        """
        Create a plot of the orbit.

        Parameters:
        marker: Marker style for the scatter plot.
        s: Size of markers for the scatter plot.
        figsize: Size of the figure.
        dpi: Dots per inch.
        fontsize: Font size for the title and labels.
        legend_loc: Location of the legend.
        legend_bbox_to_anchor: The bbox that the legend will be anchored to.
        linewidth (float): Width of the lines in the plot.

        Returns:
        fig: The created figure.
        ax: The created axes.
        """
        y_reshaped = self.y.reshape(len(self.time), -1, 4)
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        fig.suptitle(f"{self.name}", fontsize=fontsize)
        ax.set_xlabel("X", fontsize=fontsize)
        ax.set_ylabel("Y", fontsize=fontsize)

        custom_colors = ["b", "r", "k"]
        color = iter(
            plt.cm.rainbow(
                np.linspace(0, 1, max(0, y_reshaped.shape[1] - len(custom_colors)))
            )
        )

        for i in range(y_reshaped.shape[1]):
            if i < len(custom_colors):
                c_hex = custom_colors[i]
            else:
                c = next(color)
                c_hex = to_hex(c)
            ax.plot(
                y_reshaped[:, i, 0],
                y_reshaped[:, i, 2],
                color=c_hex,
                label=f"cuerpo {i+1}",
                alpha=0.5,
                linewidth=linewidth,
            )
            ax.scatter(
                y_reshaped[-1, i, 0], y_reshaped[-1, i, 2], marker=marker, s=s, c=c_hex
            )

        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height])
        ax.legend(
            loc=legend_loc, bbox_to_anchor=legend_bbox_to_anchor, fontsize=fontsize
        )
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

    def calculate_kinetic_energy(self):
        y_reshaped = self.y.reshape(len(self.time), -1, 4)
        velocities = y_reshaped[:, :, [1, 3]]
        kinetic_energy = 0.5 * self.M * np.sum(velocities**2, axis=2)
        total_kinetic_energy = np.sum(kinetic_energy, axis=1)
        return total_kinetic_energy

    def calculate_potential_energy(self):
        y_reshaped = self.y.reshape(len(self.time), -1, 4)
        positions = y_reshaped[:, :, [0, 2]]
        differences = positions[:, :, np.newaxis, :] - positions[:, np.newaxis, :, :]
        distances = np.linalg.norm(differences, axis=-1)
        pairwise_potential_energy = np.divide(
            -G * self.M[np.newaxis, :, np.newaxis] * self.M[np.newaxis, :],
            distances,
            where=distances != 0,
        )
        total_potential_energy = np.sum(pairwise_potential_energy, axis=(1, 2))
        return total_potential_energy

    def calculate_total_energy(self):
        self.kinetic_energy = self.calculate_kinetic_energy()
        self.potential_energy = self.calculate_potential_energy()
        self.total_energy = self.kinetic_energy + self.potential_energy

    def _create_plot_total_energy(self):
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

    def _create_plot_energy(self):
        fig, ax1 = plt.subplots()
        fig.suptitle("Tiempo vs Energia Total")

        ax1.plot(self.time, self.kinetic_energy, color="r", alpha=0.5)
        ax1.set_xlabel("Tiempo")
        ax1.set_ylabel("Energia Cinetica", color="r")
        ax1.tick_params(axis="y", labelcolor="r")

        ax2 = ax1.twinx()
        ax2.plot(self.time, self.potential_energy, color="b")
        ax2.set_ylabel("Energia Potencial", color="b")
        ax2.tick_params(axis="y", labelcolor="b")

        ax3 = ax1.twinx()
        # Move the third y-axis to the right by 60 points
        ax3.spines["right"].set_position(("outward", 60))
        ax3.plot(self.time, self.total_energy, color="k", alpha=0.5)
        ax3.set_ylabel("Energia Total", color="k")
        ax3.tick_params(axis="y", labelcolor="k")

        fig.tight_layout()
        return fig, (ax1, ax2, ax3)

    def plot_energy(self, route=None):
        if route is None:
            route = self.save_route_images
        fig, _ = self._create_plot_energy()
        fig.savefig(
            os.path.join(
                route, f"{self.name}_{self.ODESolver}_{self.method}_energy.png"
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
        self.plot_energy()
        self.save_error()


if "__main__" == __name__:
    init_setup = {
        "name": "test",
        "G": 1,
        "M": [1, 1, 1],
        "y1": [-0.97000436, 0.4662036850, 0.24208753, 0.4323657300],  # x1, vx1, y1, vy1
        "y2": [0.0, -0.933240737, 0.0, -0.86473146],
        "y3": [0.97000436, 0.4662036850, -0.24208753, 0.4323657300],
        "T": 6.3259,
        "h": 1e-3,
    }

    bdrk = BodySystem(init_setup=init_setup, ODESolver="rk", method="four")
    t, y = bdrk.run_simulation()
    # bdrk.save_simulation()
    bdrk.plot_energy()
    # bdrk.load_simulation()
    # print(bdrk.time)
    # print(bdrk.y)
    bdrk.plot_orbit()
    bdrk.plot_total_energy()
    # bdsym = BodySystem(ODESolver="sym", method="verlet")
    # t, y = bdsym.run_simulation()
    # bdsym.plot_orbit()
