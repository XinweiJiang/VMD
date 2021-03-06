Variational Mode Decomposition (VMD) and Its Variants 
========
The orginal VMD code:  VMD.m 

    K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans. on Signal Processing

The Multivariate Variational Mode Decomposition code:  MVVMD.m

    N. Rehman, H. Aftab, Multivariate Variational Mode Decomposition, arXiv:1907.04509, 2019.

Our works:  MVMD.p,  STMVMD.p,  MAC.p, MVMD.pyd, STMVMD.pyd. 
Only pcodes for Matlab R2016a and pydcodes for Python 3.6.5 are available now.
Please note: we only permit to use these programs to verify our paper, "Multi-dimensional Variational Mode Decomposition and Its Short-time Counterpart".
Other purposes are not permitted until further notice.
If you have any questions regarding the above codes, please contact me at liushuaishuai_hit@163.com.

Input and Parameters:
=======
signal - the time domain signal to be decomposed

alpha - the balancing parameter of the data-fidelity constraint

tau - time-step of the dual ascent ( pick 0 for noise-slack )

K - the number of modes to be recovered

DC - true if the first mode is put and kept at DC (0-freq)

init - 0 = all omegas start at 0

    - 1 = all omegas start uniformly distributed   
    
    - 2 = all omegas initialized randomly
tol - tolerance of convergence criterion; typically around 1e-6

winLen - the number of analysis points of a sliding window

overlap - the number of overlap points of adjacent windows
