import os
import time
import uuid
from collections.abc import Callable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import tomli_w
import tomllib
from matplotlib.colors import to_hex
from ODESolvers.methods import rk_methods, sym_methods
from ODESolvers.rungekutta import RKMethod
from ODESolvers.symplectic import SymIntegrator
from routes import save_route_data, save_route_images

mpl.rcParams["figure.figsize"] = (9.6, 7.2)

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

        file_name = f"{self.name}_{self.ODESolver}_{self.method}"
        self.path_data = os.path.join(self.save_route_data, file_name)

        file_name = f"{self.name}_{self.ODESolver}_{self.method}.toml"
        self.file_path_toml = os.path.join(self.save_route_data, file_name)

    def save_setup_to_toml(self):
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
        self.calculate_total_energy()
        self.calculate_linear_momentum()
        return self.time, self.y

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

    def save_simulation(self, file_format: str = "npy") -> None:
        simulation_data = np.hstack((self.time[:, np.newaxis], self.y))
        try:
            if file_format == "npy":
                with open(f"{self.path_data}.{file_format}", "wb") as file:
                    np.save(file, simulation_data)
            elif file_format == "txt":
                np.savetxt(f"{self.path_data}.{file_format}", simulation_data)
            else:
                raise ValueError("Invalid file format")
        except OSError:
            print(f"Error saving data to {self.path_data}")

    def load_simulation(self, file_format: str = "npy") -> None:
        try:
            if file_format == "npy":
                with open(f"{self.path_data}.{file_format}", "rb") as file:
                    data = np.load(file)
                    self.time = data[:, 0]
                    self.y = data[:, 1:]
            elif file_format == "txt":
                data = np.loadtxt(self.path_data)
                self.time = data[:, 0]
                self.y = data[:, 1:]
        except Exception as e:
            print(f"Error loading file: {e}")

    def calculate_kinetic_energy(self):
        y_reshaped = self.y.reshape(len(self.time), -1, 4)
        velocities = y_reshaped[:, :, [1, 3]]
        kinetic_energy = 0.5 * self.M * np.sum(velocities**2, axis=2)
        total_kinetic_energy = np.sum(kinetic_energy, axis=1)
        self._save_error_metric("Error Energia Cinetica", total_kinetic_energy)
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
        self._save_error_metric("Error Energia Potencial", total_potential_energy)
        return total_potential_energy

    def calculate_total_energy(self):
        kinetic_energy = self.calculate_kinetic_energy()
        potential_energy = self.calculate_potential_energy()
        total_energy = kinetic_energy + potential_energy
        self._save_error_metric("Error Energia Total", total_energy)
        return total_energy

    def calculate_angular_momentum(self):
        y_reshaped = self.y.reshape(len(self.time), -1, 4)
        velocities = y_reshaped[:, :, [1, 3]]
        positions = y_reshaped[:, :, [0, 2]]
        angular_momentum = np.sum(self.M * np.cross(positions, velocities), axis=1)
        self._save_error_metric("Error Energia Momento Angular", angular_momentum)
        return angular_momentum

    def calculate_linear_momentum(self):
        y_reshaped = self.y.reshape(len(self.time), -1, 4)
        velocities = y_reshaped[:, :, [1, 3]]
        linear_momentum = np.sum(self.M.reshape(3, 1) * velocities, axis=1)
        linear_momentum_x = linear_momentum[:, 0]
        linear_momentum_y = linear_momentum[:, 1]
        self._save_error_metric("Error Energia Momento Lineal X", linear_momentum_x)
        self._save_error_metric("Error Energia Momento Lineal Y", linear_momentum_y)
        return linear_momentum_x, linear_momentum_y

    def generate_orbit_figure(
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

    def generate_energy_figure(self, show_total_energy=False, alpha=1):
        kinetic_energy = self.calculate_kinetic_energy()
        potential_energy = self.calculate_potential_energy()
        total_energy = kinetic_energy + potential_energy
        fig, ax1 = plt.subplots()
        fig.suptitle("Tiempo vs Energia Total")

        ax1.plot(self.time, kinetic_energy, color="r", alpha=alpha)
        ax1.set_xlabel("Tiempo")
        ax1.set_ylabel("Energia Cinetica", color="r")
        ax1.tick_params(axis="y", labelcolor="r")

        ax2 = ax1.twinx()
        ax2.plot(self.time, potential_energy, color="b", alpha=alpha)
        ax2.set_ylabel("Energia Potencial", color="b")
        ax2.tick_params(axis="y", labelcolor="b")

        if show_total_energy:
            ax3 = ax1.twinx()
            # Move the third y-axis to the right by 60 points
            ax3.spines["right"].set_position(("outward", 60))
            ax3.plot(self.time, total_energy, color="k", alpha=alpha)
            ax3.set_ylabel("Energia Total", color="k")
            ax3.tick_params(axis="y", labelcolor="k")
            return fig, (ax1, ax2, ax3)

        fig.tight_layout()
        return fig, (ax1, ax2)

    def generate_total_energy_figure(self):
        total_energy = self.calculate_total_energy()
        fig, ax = plt.subplots()
        fig.suptitle("Tiempo vs Energia Total")
        ax.plot(self.time, total_energy)
        ax.set_xlabel("Tiempo")
        ax.set_ylabel("Energia Total")
        return fig, ax

    def generate_angular_momentum_figure(self):
        angular_momentum = self.calculate_angular_momentum()
        fig, ax = plt.subplots()
        fig.suptitle("Tiempo vs Momento Angular")
        ax.plot(self.time, angular_momentum)
        ax.set_xlabel("Tiempo")
        ax.set_ylabel("Momento angular")
        return fig, ax

    def generate_linear_momentum_figure(self):
        linear_momentum_x, linear_momentum_y = self.calculate_linear_momentum()
        fig, ax1 = plt.subplots()
        fig.suptitle("Tiempo vs Momento Lineal")

        ax1.plot(self.time, linear_momentum_x, color="r")
        ax1.set_xlabel("Tiempo")
        ax1.set_ylabel("Momemento Lineal en x", color="r")
        ax1.tick_params(axis="y", labelcolor="r")

        ax2 = ax1.twinx()
        ax2.plot(self.time, linear_momentum_y, color="b")
        ax2.set_xlabel("Tiempo")
        ax2.set_ylabel("Momemento Lineal en y", color="b")
        ax2.tick_params(axis="y", labelcolor="b")

        fig.tight_layout()
        return fig, (ax1, ax2)

    def save_orbit_figure(self, route=None):
        """
        Plot the orbit and save the figure.

        Parameters:
        route (str): The path where the figure will be saved. If None, the default path is used.

        Returns:
        None
        """
        if route is None:
            route = self.save_route_images

        if not os.path.exists(route):
            os.makedirs(route, exist_ok=True)

        fig, _ = self.generate_orbit_figure()
        filename = f"{self.name}_{self.ODESolver}_{self.method}.png"
        filepath = os.path.join(route, filename)
        fig.savefig(filepath)
        plt.close(fig)

    def save_energy_figure(self, route=None):
        if route is None:
            route = self.save_route_images

        if not os.path.exists(route):
            os.makedirs(route, exist_ok=True)

        fig, _ = self.generate_energy_figure()
        filename = f"{self.name}_{self.ODESolver}_{self.method}_energy.png"
        filepath = os.path.join(route, filename)
        fig.savefig(filepath)
        plt.close(fig)

        fig, _ = self.generate_energy_figure(True, alpha=0.5)
        filename = f"{self.name}_{self.ODESolver}_{self.method}_energy_t.png"
        filepath = os.path.join(route, filename)
        fig.savefig(filepath)
        plt.close(fig)

    def save_total_energy_figure(self, route=None):
        if route is None:
            route = self.save_route_images

        if not os.path.exists(route):
            os.makedirs(route, exist_ok=True)

        fig, _ = self.generate_total_energy_figure()
        filename = f"{self.name}_{self.ODESolver}_{self.method}_total_energy.png"
        filepath = os.path.join(route, filename)
        fig.savefig(filepath)
        plt.close(fig)

    def save_angular_momentum_figure(self, route=None):
        if route is None:
            route = self.save_route_images

        if not os.path.exists(route):
            os.makedirs(route, exist_ok=True)

        fig, _ = self.generate_angular_momentum_figure()
        filename = f"{self.name}_{self.ODESolver}_{self.method}_angular_momentum.png"
        filepath = os.path.join(route, filename)
        fig.savefig(filepath)
        plt.close(fig)

    def save_linear_momentum_figure(self, route=None):
        if route is None:
            route = self.save_route_images

        if not os.path.exists(route):
            os.makedirs(route, exist_ok=True)

        fig, _ = self.generate_linear_momentum_figure()
        filename = f"{self.name}_{self.ODESolver}_{self.method}_linear_momentum.png"
        filepath = os.path.join(route, filename)
        fig.savefig(filepath)
        plt.close(fig)

    # TODO: Revisar el error abosluto y relativo
    def _error(self, data):
        reference_value = data[0]
        if reference_value == 0:
            reference_value = EPS
        absolute_errors = np.abs(data[1:] - reference_value)
        relative_errors = absolute_errors / np.abs(reference_value)
        return np.mean(absolute_errors), np.mean(relative_errors)

    def _save_error_metric(self, metric_name, value):
        with open(self.file_path_toml, "rb") as f:
            data = tomllib.load(f)

        if "metrics" not in data:
            data["metrics"] = {}

        data["metrics"][metric_name] = self._error(value)

        with open(self.file_path_toml, "wb") as f:
            tomli_w.dump(data, f)

    def plot(self):
        self.save_orbit_figure()
        self.save_energy_figure()
        self.save_total_energy_figure()
        self.save_angular_momentum_figure()
        self.save_linear_momentum_figure()


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
    bdrk.load_simulation()
    # t, y = bdrk.run_simulation()
    bdrk.plot()
