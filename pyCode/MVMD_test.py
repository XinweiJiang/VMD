# Copyright(c) 2020 Shuaishuai Liu.All rights reserved.

# We only permit to use these programs to verify our
# paper, "Multi-dimensional Variational Mode Decomposition and Its Short-time Counterpart".
# Other purposes are not permitted until further notice.

import numpy as np
import matplotlib.pyplot as plt
from mvmd import mvmd


fs = 1024  # sample frequency
t = np.linspace(0, 1, fs)  # discrete time step
fc = np.array([[20], [124], [186]])  # frequency contents
s = np.sin(2 * np.pi * fc * t)  # IMF

rd = np.random.RandomState(1)  # random seed

B = rd.normal(0, 1, [2, 3])  # mixing matrix
x = 10*np.matrix(B)*np.matrix(s)  # mixing signals

# solution parameters
alpha = 2000
tau = 0
K = 3
DC = 0
init = 1
tol = 1e-6

u, omega, B = mvmd(x, alpha, tau, K, DC, init, tol)


# plot the iteration history of the center frequency
plt.figure(1)
Index = np.linspace(1, omega.shape[0], omega.shape[0])
plt.plot(omega[:, 0]*fs, Index, color='red')
plt.plot(omega[:, 1]*fs, Index, color='green')
plt.plot(omega[:, 2]*fs, Index, color='blue')
plt.title('The iteration history of the center frequency')
plt.xlabel('frequency / Hz')
plt.ylabel('iteration step')
plt.legend(('$f_1$', '$f_2$', '$f_3$'), loc='upper right')
plt.grid(True)

# plt.show()



# --------------------------------
# input data x(t)
# plot x1(t)
plt.figure(2)
plt.subplot(421)
plt.plot(t, x[0, :].T, '-k')
plt.title('$x_1$(t)')
plt.xticks([])

# plot x2(t)
plt.subplot(422)
plt.plot(t, x[1, :].T, '-k')
plt.title('$x_2$(t)')
plt.xticks([])

# --------------------------------
# decomposition results of channel x1(t)
# plot u1(t) decomposed by x1(t)
plt.subplot(423)
plt.plot(t, B[0, 0]*u[0, :], '-k')
plt.ylabel('$u_1$(t)')
plt.xticks([])

# plot u2(t) decomposed by x1(t)
plt.subplot(425)
plt.plot(t, B[0, 1]*u[1, :], '-k')
plt.ylabel('$u_2$(t)')
plt.xticks([])

# plot u1(t) decomposed by x1(t)
plt.subplot(427)
plt.plot(t, B[0, 2]*u[2, :], '-k')
plt.ylabel('$u_3$(t)')
plt.xlabel('time / Sec')


# --------------------------------
# decomposition results of channel  x2(t)
# plot u1(t) decomposed by x2(t)
plt.subplot(424)
plt.plot(t, B[1, 0]*u[0, :], '-k')
plt.xticks([])

# plot u2(t) decomposed by x2(t)
plt.subplot(426)
plt.plot(t, B[1, 1]*u[1, :], '-k')
plt.xticks([])

# plot u1(t) decomposed by x2(t)
plt.subplot(428)
plt.plot(t, B[1, 2]*u[2, :], '-k')
plt.xlabel('time / Sec')
plt.grid(True)

plt.show()