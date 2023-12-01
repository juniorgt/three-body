import os

import numpy as np
from routes import save_route_data


def txt_to_npy(name: str, ODESolver: str, method: str) -> None:
    dir_save_data = os.path.join(save_route_data, name, f"{ODESolver}_{method}")
    file_name = f"{name}_{ODESolver}_{method}"
    path_data = os.path.join(dir_save_data, file_name)

    try:
        data = np.loadtxt(path_data)
    except Exception as e:
        print(f"Error loading file: {e}")

    try:
        with open(f"{path_data}.npy", "wb") as file:
            np.save(file, data)
    except Exception as e:
        print(f"Error loading file: {e}")
