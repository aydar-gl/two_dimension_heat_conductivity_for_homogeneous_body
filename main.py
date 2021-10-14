import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Nx = 100
Ny = 100
t_end = 60
L = 0.5
H = 0.5
lamda = 46
ro = 7800
c = 460
Th = 80
Tc = 30
T0 = 5

T = np.zeros((Nx, Ny), dtype="double")
alfa = np.zeros(Nx, dtype="double")
beta = np.zeros(Nx, dtype="double")

hx = L/(Nx-1)
hy = H/(Ny-1)
a = lamda/(ro*c)
tau = t_end/100.0
T[:,:] = T0
timer = 0

while timer < t_end:
    timer += tau
    for j in range(0, Ny):
        alfa[0] = 0.0
        beta[0] = Th
        for i in range(1, Nx-1):
            ai = lamda/np.square(hx)
            bi = 2.0*lamda/np.square(hx)+ro*c/tau
            ci = lamda/np.square(hx)
            fi = -ro*c*T[i,j]/tau
            alfa[i] = ai / (bi - ci * alfa[i - 1])
            beta[i] = (ci * beta[i - 1] - fi) / (bi - ci * alfa[i - 1])
        T[Nx-1, j] = Tc
        for i in reversed(range(0, Nx-1)):
            T[i, j] = alfa[i] * T[i + 1, j] + beta[i]
    for i in range(1, Nx-1):
        alfa[0] = 2.0 * a * tau / (2.0 * a * tau + np.square(hy))
        beta[0] = np.square(hy) * T[i, 1] / (2.0 * a * tau + np.square(hy))
        for j in range(1, Ny-1):
            ai = lamda / np.square(hy)
            bi = 2.0 * lamda / np.square(hy) + ro * c / tau
            ci = lamda / np.square(hy)
            fi = -ro * c * T[i, j] / tau
            alfa[j] = ai / (bi - ci * alfa[j - 1])
            beta[j] = (ci * beta[j - 1] - fi) / (bi - ci * alfa[j - 1])
        T[i, Ny-1] = (2.0 * a * tau * beta[Ny - 2] + np.square(hy) * T[i, Ny-1]) / (2.0 * a * tau * (1.0 - alfa[Ny - 2]) + np.square(hy))
        for j in reversed(range(0, Ny-1)):
            T[i, j] = alfa[j] * T[i, j + 1] + beta[j]

f = lambda x, y: T[x,y]
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection = '3d')
X1 = np.arange(0, Nx, 1)
Y1 = np.arange(0, Ny, 1)
x1, y1 = np.meshgrid(X1, Y1)
z1 = f(x1, y1)
surf = ax.plot_surface(x1,y1,z1, rstride = 1, cstride = 1, cmap = cm.viridis)
plt.show()
