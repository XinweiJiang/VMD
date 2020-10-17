# Copyright(c) 2020 Shuaishuai Liu.All rights reserved.

# We only permit to use these programs to verify our
# paper, "Multi-dimensional Variational Mode Decomposition and Its Short-time Counterpart".
# Other purposes are not permitted until further notice.

import numpy as np
import matplotlib.pyplot as plt
from STMVMD import stmvmd
from scipy.signal import hilbert
from tftb.processing import inst_freq

fs = 1024  # sample frequency
t = np.linspace(0, 10, 10*fs)  # discrete time step
fc = np.array([[50], [124], [186]])  # frequency contents
s = np.sin(2 * np.pi * (fc * t + 0.5*np.square(t) + 20/np.pi*np.sin(0.9*np.pi*t)))  # IMF

rd = np.random.RandomState(1)  # random seed

B = rd.normal(0, 1, [2, 3])  # mixing matrix
x = 10*np.matrix(B)*np.matrix(s)  # mixing signals

# solution parameters
alpha = 1000
tau = 0
K = 3
DC = 0
init = 1
tol = 1e-6
winLen = 512
overlap = 512-128

u, T, omega, Tc, B = stmvmd(x, alpha, tau, K, DC, init, tol, winLen, overlap)

print(u.shape)

# --------------------------------
# plot the iteration history of the center frequency
plt.figure(1)
plt.plot(t[Tc], omega[:, 0]*fs, color='red')
plt.plot(t[Tc], omega[:, 1]*fs, color='green')
plt.plot(t[Tc], omega[:, 2]*fs, color='blue')
plt.title('The time-frequency curve of the center frequency \n Obtained by plotting (center time, center frequency) of all windows')
plt.ylabel('frequency / Hz')
plt.xlabel('time / Sec')
plt.legend(('$u_1$', '$u_2$', '$u_3$'), loc='upper right')
plt.grid(True)


# --------------------------------
# decomposition results
# plot u1(t)
plt.figure(2)
plt.subplot(331)
u1 = u[0, :]
plt.plot(t, u1, '-k')
plt.ylabel('$u_1$(t)')
plt.title('IMF $u_k$')
plt.xticks([])

# plot the instantaneous amplitude of u1
plt.subplot(332)
h = hilbert(u1)
A = np.abs(h)
plt.plot(t, A, '-k')
plt.title('Amplitude of $u_k$')
plt.xticks([])

# plot the instantaneous frequency of u1
plt.subplot(333)
instf, timestamps = inst_freq(h)
plt.plot(timestamps/fs, instf*fs, '.k')
plt.title('Frequency of $u_k$')
plt.xticks([])


# plot u2(t)
plt.subplot(334)
u2 = u[1, :]
plt.plot(t, u2, '-k')
plt.ylabel('$u_2$(t)')
plt.xticks([])

# plot the instantaneous amplitude of u2
plt.subplot(335)
h = hilbert(u2)
A = np.abs(h)
plt.plot(t, A, '-k')
plt.xticks([])

# plot the instantaneous frequency of u2
plt.subplot(336)
instf, timestamps = inst_freq(h)
plt.plot(timestamps/fs, instf*fs, '.k')
plt.xticks([])

# plot u3(t)
plt.subplot(337)
u3 = u[2, :]
plt.plot(t, u3, '-k')
plt.ylabel('$u_3$(t)')
plt.xlabel('time / Sec')

# plot the instantaneous amplitude of u3
plt.subplot(338)
h = hilbert(u3)
A = np.abs(h)
plt.plot(t, A, '-k')
plt.xlabel('time / Sec')

# plot the instantaneous frequency of u3
plt.subplot(339)
instf, timestamps = inst_freq(h)
plt.plot(timestamps/fs, instf*fs, '.k')
plt.xlabel('time / Sec')


# plt.grid(True)

plt.show()
