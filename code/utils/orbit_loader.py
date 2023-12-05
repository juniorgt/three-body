import tomllib


def load_toml(file_path: str) -> dict:
    try:
        with open(file_path, "rb") as file:
            data = tomllib.load(file)
    except FileNotFoundError:
        print(f"No se encontró el archivo en la ruta {file_path}")

    return data


def get_orbit(name_orbit: str, file: str = "orbits.toml") -> dict:
    orbits = load_toml(file)
    if name_orbit not in orbits:
        raise ValueError(
            f"La órbita '{name_orbit}' no se encontró en el archivo de órbitas"
        )
    return orbits[name_orbit]


def get_all_orbits(file: str = "orbits.toml") -> dict:
    orbits = load_toml(file)
    return orbits
