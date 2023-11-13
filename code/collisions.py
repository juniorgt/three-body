import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Function to calculate delta_r
def delta_r(
    coords, masses, nBodies, G
):  # esta es una versión generalizada de lo que teníamos antes para nCuerpos, diferentes masas y diferentes valores de G
    x, y, vx, vy = coords.copy()  # esto crea una copia física en memoria
    delta = coords.copy()

    for n in range(nBodies):
        xn, yn = (
            x[n],
            y[n],
        )  # asignarles nuevas variables aquí por eficiencia computacional, minimizar las llamadas de acceso a memoria
        delta_vx, delta_vy = 0.0, 0.0
        for i in range(
            nBodies
        ):  # generalizando para el problema posterior de n-cuerpos
            if i != n:  # sólo calcular si no es mismo
                sep = np.sqrt(
                    (xn - x[i]) ** 2 + (yn - y[i]) ** 2
                )  # Distancia euclidiana
                delta_vx -= (
                    G * masses[i] * (xn - x[i]) / sep**3
                )  # cambio de velocidad de cada masa sobre la masa n
                delta_vy -= G * masses[i] * (yn - y[i]) / sep**3
        delta[2, n] = delta_vx  # cambio de velocidad = a*t
        delta[3, n] = delta_vy
    delta[0] = vx  # cambio de posición = v*t
    delta[1] = vy
    return delta


# Function for one RK4 step
def step(
    coords, masses, delta_t, nBodies=3, G=6.67408313131313e-11
):  # 1 paso RK4 para las coordenadas de cada cuerpo, muta las coordenadas
    k1 = delta_t * delta_r(coords, masses, nBodies, G)
    k2 = delta_t * delta_r(coords + k1 / 2, masses, nBodies, G)
    k3 = delta_t * delta_r(coords + k2 / 2, masses, nBodies, G)
    k4 = delta_t * delta_r(coords + k3, masses, nBodies, G)
    coords += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6
    return coords


# Initial conditions setup
# #Condiciones iniciales -- ¡puedes modificarlas para hacer las animaciones de tres cuerpos que quieras!
M = [1, 1, 1]  # m1=m2=m3
x = [-0.97000436, 0.0, 0.97000436]  # x1 = -x3, x2 = 0
y = [0.24208753, 0.0, -0.24208753]  # y1 = -y3, y2 = 0
vx = [0.4662036850, -0.933240737, 0.4662036850]  # v1x = v3x
vy = [0.4323657300, -0.86473146, 0.4323657300]  # v1y = v3y
Ei = (
    -1 / np.sqrt((2 * 0.97000436) ** 2 + (2 * 0.24208753) ** 2)
    - 2 / np.sqrt(0.97000436**2 + 0.24208753**2)
    + 0.5 * sum(np.array(vx) ** 2 + np.array(vy) ** 2)
)  # r23 = r12
coords = np.array([x, y, vx, vy])
time = np.linspace(0, 6.3259, 1001)
Δt = time[1] - time[0]

### coords
# [[x1, x2, x3],
#  [y1, y2, y3],
#  [vx1, vx2, vx3],
#  [vy1, vy2, vy3]]

X = np.zeros((3, len(time)))
Y = np.zeros((3, len(time)))
VX = np.zeros((3, len(time)))
VY = np.zeros((3, len(time)))

for i in range(len(time)):
    coords = step(coords, M, Δt, 3, 1)  # evolucionar el sistema en cada paso temporal
    X[:, i] = coords[0]  # mantener la posición en cada paso temporal
    Y[:, i] = coords[1]
    VX[:, i] = coords[2]
    VY[:, i] = coords[3]


def detect_collisions_escape(coords, delta_t, max_sep):
    x, y, vx, vy = coords
    V = np.sqrt(np.square(vx) + np.square(vy))
    R = V * delta_t
    collision = False
    collision_inds = None
    escape = False
    escape_ind = None

    for n in range(len(x)):
        rn = R[n]
        xn = x[n]
        yn = y[n]

        for i in range(len(x)):
            if i != n:
                min_sep = rn + R[i]
                sep = np.sqrt((xn - x[i]) ** 2 + (yn - y[i]) ** 2)

                if sep < min_sep:
                    collision = True
                    collision_inds = (n, i)
                    return collision, collision_inds, escape, escape_ind
                elif sep > max_sep:
                    escape = True
                    escape_ind = n
                    return collision, collision_inds, escape, escape_ind

    return collision, collision_inds, escape, escape_ind


def perturb_test(delta_max=0.1):
    # Initial conditions setup
    M = [1, 1, 1]  # m1 = m2 = m3
    x = [-0.97000436, 0.0 + random.uniform(-1, 1) * np.sqrt(delta_max), 0.97000436]
    y = [0.24208753, 0.0 + random.uniform(-1, 1) * np.sqrt(delta_max), -0.24208753]
    vx = [0.4662036850, -0.933240737, 0.4662036850]
    vy = [0.4323657300, -0.86473146, 0.4323657300]

    coords = np.array([x, y, vx, vy])  # array of initial conditions
    init_coords = np.copy([x[1], y[1]])  # keep track for output later
    time = np.linspace(0, 6.3259 * 10, num=1001)
    delta_t = (
        time[1] - time[0]
    )  # the periodicity of this system is ~6.3259 in these units
    max_sep = 10 * np.sqrt(
        (2 * 0.97000436) ** 2 + (2 * 0.24208753) ** 2
    )  # biggest distance apart is at start between bodies 1-3
    collision, collision_inds, escape, escape_ind = False, None, False, None

    for i in range(len(time)):
        # Implement a step! function that updates the coordinates
        coords = step(coords, M, delta_t, 3, 1)

        # Implement a detectCollisionsEscape function to check for collisions or escapes
        collision, collision_inds, escape, escape_ind = detect_collisions_escape(
            coords, delta_t, max_sep
        )

        if collision or escape:
            return collision, collision_inds, escape, escape_ind, init_coords, time[i]

    return collision, collision_inds, escape, escape_ind, init_coords, time[-1]


N = 10000  # number of random trials, make this number higher to get a prettier/more accurate plot (but it also takes longer to run...)
column_names = ["collision", "escape", "T", "x", "y", "r"]
column_types = [bool, bool, float, float, float, float]

data = {name: [None] * N for name in column_names}
df = pd.DataFrame(data)

for i in range(1, N + 1):
    if i % (N // 1000) == 0:
        print(f"{i / N * 100:.2f} % complete", end="\r")

    collision, collision_inds, escape, escape_ind, init_coords, T = perturb_test()
    df.at[i - 1, "collision"] = collision
    df.at[i - 1, "escape"] = escape
    df.at[i - 1, "T"] = T
    df.at[i - 1, "x"] = init_coords[0]
    df.at[i - 1, "y"] = init_coords[1]
    df.at[i - 1, "r"] = np.sqrt(init_coords[0] ** 2 + init_coords[1] ** 2)

print()

# Visualize results
# variables x, y -> Condiciones iniciales
plt.figure(figsize=(12, 4))
plt.scatter(df["x"], df["y"], c=df["T"], s=2, alpha=0.8)
plt.plot(
    X[1, :],
    Y[1, :],
    linestyle="--",
    color="black",
    lw=2,
    alpha=0.7,
    label="órbita ∞ estable: estrellas + flechas = condiciones iniciales estables",
)
plt.scatter(x, y, marker="*", s=200, edgecolor="none", c=["r", "m", "b"])

quiver_scale = 0.2
plt.quiver(X[0, 0], Y[0, 0], VX[0, 0], VY[0, 0], minshaft=3, color="r")
plt.quiver(X[1, 0], Y[1, 0], VX[1, 0], VY[1, 0], minshaft=3, color="m")
plt.quiver(X[2, 0], Y[2, 0], VX[2, 0], VY[2, 0], minshaft=3, color="b")

plt.title("Duración de la simulación en función de la posición inicial del cuerpo 2")
plt.xlabel("x")
plt.ylabel("y")
# plt.gca().set_aspect('equal')
plt.legend()
# plt.colorbar()

# Save the plot as a PNG
plt.savefig("chaos.png")
plt.close()

# Display the saved image
image = plt.imread("chaos.png")
plt.figure(figsize=(image.shape[1] / 100, image.shape[0] / 100))
plt.imshow(image)
plt.axis("off")
plt.show()
