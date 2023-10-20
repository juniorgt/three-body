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


def runge_kutta_fehlberg(f, coords, masses, h, nBodies, G):
    k1 = h * f(coords, masses, nBodies, G)
    k2 = h * f(coords + (1 / 4) * k1, masses, nBodies, G)
    k3 = h * f(coords + (3 / 32) * k1 + (9 / 32) * k2, masses, nBodies, G)
    k4 = h * f(
        coords + (1932 / 2197) * k1 - (7200 / 2197) * k2 + (7296 / 2197) * k3,
        masses,
        nBodies,
        G,
    )
    k5 = h * f(
        coords + (439 / 216) * k1 - 8 * k2 + (3680 / 513) * k3 - (845 / 4104) * k4,
        masses,
        nBodies,
        G,
    )
    k6 = h * f(
        coords
        - (8 / 27) * k1
        + 2 * k2
        - (3544 / 2565) * k3
        + (1859 / 4104) * k4
        - (11 / 40) * k5,
        masses,
        nBodies,
        G,
    )

    coords += (
        (16 / 135) * k1
        + (6656 / 12825) * k3
        + (28561 / 56430) * k4
        - (9 / 50) * k5
        + (2 / 55) * k6
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
