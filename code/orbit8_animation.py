import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation

def delta_r(coords, masses, nBodies, G):  
    x, y, vx, vy = coords.copy()
    delta = coords.copy()
    for n in range(nBodies):
        xn, yn = x[n], y[n]
        delta_vx, delta_vy = 0., 0.
        for i in range(nBodies): 
            if i != n:
                sep = np.sqrt((xn - x[i]) ** 2 + (yn - y[i]) ** 2)
                delta_vx -= G * masses[i] * (xn - x[i]) / sep ** 3
                delta_vy -= G * masses[i] * (yn - y[i]) / sep ** 3
        delta[2, n] = delta_vx 
        delta[3, n] = delta_vy
    delta[0] = vx
    delta[1] = vy
    return delta


def step(coords, masses, delta_t, nBodies=3, G=6.67408313131313e-11): 
    k1 = delta_t * delta_r(coords, masses, nBodies, G)
    k2 = delta_t * delta_r(coords + k1 / 2, masses, nBodies, G)
    k3 = delta_t * delta_r(coords + k2 / 2, masses, nBodies, G)
    k4 = delta_t * delta_r(coords + k3, masses, nBodies, G)
    coords += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6
    return coords

M = [1, 1, 1] #m1=m2=m3
x = [-0.97000436, 0., 0.97000436]  #x1 = -x3, x2 = 0
y = [0.24208753, 0., -0.24208753] #y1 = -y3, y2 = 0
vx = [0.4662036850, -0.933240737, 0.4662036850] #v1x = v3x
vy = [0.4323657300, -0.86473146, 0.4323657300] #v1y = v3y
Ei = -1 / np.sqrt((2 * 0.97000436) ** 2 + (2 * 0.24208753) ** 2) - 2 / np.sqrt(0.97000436 ** 2 + 0.24208753 ** 2) + 0.5 * sum(np.array(vx) ** 2 + np.array(vy) ** 2) #r23 = r12
coords = np.array([x, y, vx, vy])
time = np.linspace(0, 6.3259, 1000)
delta_t = time[1] - time[0]


X = np.zeros((3, len(time)))
Y = np.zeros((3, len(time)))
VX = np.zeros((3, len(time)))
VY = np.zeros((3, len(time)))

for i in range(len(time)):
    coords = step(coords, M, delta_t, 3, 1) 
    X[:, i] = coords[0]
    Y[:, i] = coords[1]
    VX[:, i] = coords[2]
    VY[:, i] = coords[3]



fig = plt.figure(figsize=(12, 4))
gs = GridSpec(1, 2, width_ratios=[5, 1])
ax = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]

ax[0].set_xlim(-1.1, 1.1)
ax[0].set_ylim(-0.4, 0.4)
ax[0].set_xlabel("x")
ax[0].set_ylabel("y")
ax[0].set_title("Problema de los 3 cuerpos: t = ")

bar_ax = ax[1]
bar_ax.set_ylim(0, 3e-10)
bar_ax.set_xticks([])
bar_ax.set_yticks([0, 1e-10, 2e-10, 3e-10])
bar_ax.set_yticklabels([0, '1e-10', '2e-10', '3e-10'])
bar_ax.set_ylabel("Energia perdida [%]")

line1, = ax[0].plot([], [], label="", lw=2, color='blue')
line2, = ax[0].plot([], [], label="", lw=2, color='green')
line3, = ax[0].plot([], [], label="", lw=2, color='red')
star1, = ax[0].plot([], [], label="", marker='o', color='blue', markersize=8, markeredgewidth=0)
star2, = ax[0].plot([], [], label="", marker='o', color='green', markersize=8, markeredgewidth=0)
star3, = ax[0].plot([], [], label="", marker='o', color='red', markersize=8, markeredgewidth=0)
bar = bar_ax.bar([0], 0, color='crimson', width=0.1)

def update(frame):
    ax[0].set_title("Problema de los 3 cuerpos: t = {:.3f}".format(time[frame]))
    line1.set_data(X[0, :frame], Y[0, :frame])
    line2.set_data(X[1, :frame], Y[1, :frame])
    line3.set_data(X[2, :frame], Y[2, :frame])
    star1.set_data([X[0, frame]], [Y[0, frame]])
    star2.set_data([X[1, frame]], [Y[1, frame]])
    star3.set_data([X[2, frame]], [Y[2, frame]])
    r12 = np.sqrt((X[0, frame] - X[1, frame]) ** 2 + (Y[0, frame] - Y[1, frame]) ** 2)
    r13 = np.sqrt((X[0, frame] - X[2, frame]) ** 2 + (Y[0, frame] - Y[2, frame]) ** 2)
    r23 = np.sqrt((X[2, frame] - X[1, frame]) ** 2 + (Y[2, frame] - Y[1, frame]) ** 2)
    U = -1 / r12 - 1 / r13 - 1 / r23
    K = 0.5 * np.sum(VX[:, frame] ** 2 + VY[:, frame] ** 2)
    delta_E = abs((U + K - Ei) / Ei)
    bar[0].set_height(delta_E)

# Create an animation
ani = FuncAnimation(fig, update, frames=len(time), blit=False, interval=1)
# ani.save("0test.gif", dpi=300, writer=PillowWriter(fps=25))
#ani.save("01test.mp4", writer='ffmpeg', fps=30)
plt.show()
