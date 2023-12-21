import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tomllib
from routes import save_route_data

cmap = plt.get_cmap("nipy_spectral")


def new_name(name):
    return " ".join(name.split("_"))


def buscar_carpetas(nombre: str) -> list:
    carpetas_encontradas = []
    for root, dirs, _ in os.walk(save_route_data):
        for dir in dirs:
            if nombre in dir:
                carpetas_encontradas.append(os.path.join(root, dir))
    return carpetas_encontradas


def print_latex(name_orbits: str):
    d = buscar_carpetas("Figure-8")

    data_dic = {}
    for c in d:
        for root, _, files in os.walk(c):
            if files:
                metric = None
                name_methdod = os.path.basename(root)
                path_toml = os.path.join(root, files[1])
                # print(path_toml)
                with open(path_toml, "rb") as f:
                    tm = tomllib.load(f)
                    metric = tm["metrics"]["running_time"]
                    # metric = tm["metrics"]["Error Energia Total"][0]
                if data_dic.get(name_methdod, None) is None:
                    data_dic[name_methdod] = []
                data_dic[name_methdod].append(metric)

    df = pd.DataFrame(data_dic).transpose()
    df.columns = ["1e-4", "1e-3", "1e-2", "1e-1", "1e-5"]
    df = df[["1e-1", "1e-2", "1e-3", "1e-4", "1e-5"]]
    df.index = df.index.map(new_name)
    print(df.to_latex())


def make_figure(name_orbits: str):
    d = buscar_carpetas(name_orbits)

    data_dic = {}
    for c in d:
        for root, _, files in os.walk(c):
            if files:
                metric = None
                name_methdod = os.path.basename(root)
                path_toml = os.path.join(root, files[1])
                # print(path_toml)
                with open(path_toml, "rb") as f:
                    tm = tomllib.load(f)
                    metric = tm["metrics"]["running_time"]
                if data_dic.get(name_methdod, None) is None:
                    data_dic[name_methdod] = []
                data_dic[name_methdod].append(metric)

    df = pd.DataFrame(data_dic).transpose()
    df.columns = ["1e-4", "1e-3", "1e-2", "1e-1", "1e-5"]
    df = df[["1e-1", "1e-2", "1e-3", "1e-4", "1e-5"]]
    df.index = df.index.map(new_name)
    df.to_csv("data.csv")
    # print(df.to_latex())
    df_t = df.T

    df_t.index = df_t.index.astype(float)
    df_t.sort_index(ascending=False, inplace=True)

    colores = [cmap(i) for i in np.linspace(0, 1, 24)]

    fig, ax = plt.subplots(figsize=(20, 10))
    for i, column in enumerate(df_t.columns):
        ax.plot(df_t.index, df_t[column], marker="o", label=column, color=colores[i])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_title(f"Tiempo de ejecución vs. Tamaño de paso de tiempo ({name_orbits})")
    ax.set_xlabel("Tamaño de paso")
    ax.set_ylabel("Tiempo de ejecución (segundos)")

    ax.legend()

    ax.grid(True, which="both", ls="--", linewidth=0.5)

    fig.savefig(f"results_energy/{name_orbits}_result.png")


def make_figure_energy(name_orbits: str):
    d = buscar_carpetas(name_orbits)

    data_dic = {}
    for c in d:
        for root, _, files in os.walk(c):
            if files:
                metric = None
                name_methdod = os.path.basename(root)
                path_toml = os.path.join(root, files[1])
                # print(path_toml)
                with open(path_toml, "rb") as f:
                    tm = tomllib.load(f)
                    metric = tm["metrics"]["Error Energia Total"][0]
                if data_dic.get(name_methdod, None) is None:
                    data_dic[name_methdod] = []
                data_dic[name_methdod].append(metric)

    df = pd.DataFrame(data_dic).transpose()
    df.columns = ["1e-4", "1e-3", "1e-2", "1e-1", "1e-5"]
    df = df[["1e-1", "1e-2", "1e-3", "1e-4", "1e-5"]]
    df.index = df.index.map(new_name)
    df.to_csv("data.csv")
    print(df.to_latex())
    print(df)

    pasos_tiempo = [float(col) for col in df.columns]
    pasos_tiempo.sort()
    colores = [cmap(i) for i in np.linspace(0, 1, 24)]
    # Crear la figura y los ejes
    fig, ax = plt.subplots(figsize=(20, 10))

    for i, metodo in enumerate(df.index):
        ax.plot(
            pasos_tiempo, df.loc[metodo, :], marker="o", label=metodo, color=colores[i]
        )

    # Añadir detalles al gráfico
    ax.set_title(f"Error en la Energía vs Paso de Tiempo ({name_orbits})", fontsize=14)
    ax.set_xlabel("Paso de Tiempo", fontsize=12)
    ax.set_ylabel("Error en la Energía", fontsize=12)
    ax.legend(title="Método Numérico", fontsize=10)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.tick_params(axis="both", which="major", labelsize=10)

    ax.set_xscale("log")

    ax.invert_xaxis()
    plt.savefig("grafico_tesis.png", format="png", dpi=300)

    plt.show()


if __name__ == "__main__":
    orbit_names = [
        "Figure-8",
        "Butterfly-I",
        "Butterfly-II",
        "Bumblebee",
        "Moth-I",
        "Moth-II",
        "Butterfly-III",
        "Goggles",
        "Butterfly-IV",
        "Dragonfly",
        "Yarn",
        "2a-Yin-yang-I",
        "2b-Yin-yang-I",
        "3a-Yin-yang-II",
        "3b-Yin-yang-II",
    ]
    make_figure_energy("Goggles")
    # for name in orbit_names:
    #    make_figure(name, flag_time=False)
