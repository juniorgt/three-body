def euler(f, coords, masses, h, nBodies, G):
    delta = f(coords, masses, nBodies, G)
    coords += delta * h
    return coords


def runge_kutta_4th_order(f, coords, masses, h, nBodies, G):
    k1 = h * f(coords, masses, nBodies, G)
    k2 = h * f(coords + k1 / 2, masses, nBodies, G)
    k3 = h * f(coords + k2 / 2, masses, nBodies, G)
    k4 = h * f(coords + k3, masses, nBodies, G)
    coords += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6
    return coords


def runge_kutta_5th_order(f, coords, masses, h, nBodies, G):
    k1 = h * f(coords, masses, nBodies, G)
    k2 = h * f(coords + (1 / 4) * k1, masses, nBodies, G)
    k3 = h * f(coords + (1 / 8) * k1 + (1 / 8) * k2, masses, nBodies, G)
    k4 = h * f(coords - (1 / 2) * k2 + k3, masses, nBodies, G)
    k5 = h * f(
        coords + (3 / 16) * k1 - (3 / 8) * k2 + (3 / 8) * k3 + (9 / 16) * k4,
        masses,
        nBodies,
        G,
    )
    k6 = h * f(
        coords
        - (3 / 7) * k1
        + (2 / 7) * k2
        + (12 / 7) * k3
        - (12 / 7) * k4
        + (8 / 7) * k5,
        masses,
        nBodies,
        G,
    )

    coords += (
        (7 / 90) * k1 + (32 / 90) * k3 + (12 / 90) * k4 + (32 / 90) * k5 + (7 / 90) * k6
    )

    return coords


def verlet(f, coords, masses, h, nBodies, G):
    a = f(coords, masses, nBodies, G)
    coords[0:2] += coords[2:4] * h + 0.5 * a[2:4] * h**2
    a_new = f(coords, masses, nBodies, G)
    coords[2:4] += 0.5 * (a[2:4] + a_new[2:4]) * h
    return coords


def leapfrog(f, coords, masses, h, nBodies, G):
    coords[2:4] += 0.5 * h * f(coords, masses, nBodies, G)[2:4]
    coords[0:2] += h * coords[2:4]
    coords[2:4] += 0.5 * h * f(coords, masses, nBodies, G)[2:4]
    return coords
