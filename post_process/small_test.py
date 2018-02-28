#!/usr/bin/python

""" Small FFT test
    Want to see what the FFT looks like of a simple function in a giant box.
"""

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import ticker
import config    

plt.rcParams.update(config.pars)
np.set_printoptions(threshold=np.inf)

def plot_energy_spectra():

   # Create simple function
   Nx = 2048
   n = 2.0
   xl = -n * np.pi
   xr = (Nx - 2.0) * n * np.pi / Nx
   x = np.linspace(xl, xr, Nx)
   #f = np.sin(2.0 * x)
   f = np.exp(-x * x)

   fig, ax = plt.subplots(1,1, figsize=(10,0.75*10))
   ax.plot(x,f)

   # Take FFT
   fk = np.fft.fft(f) / Nx

   # Set up wavenumbers
   kx = np.zeros(Nx)
   for ii in range(Nx//2):
      kx[ii] = float(ii)
   for ii in range(Nx//2, Nx):
      kx[ii] = float(ii-Nx)

   # Compute energy in Fourier space
   NF = Nx
   kmax = NF // 2
   Efk = np.zeros(int(kmax), dtype=np.complex_)
   k  = np.arange(0.0, kmax)
   for ii in range(Nx):
      kmag = abs(kx[ii])
      if (kmag < kmax):
         Efk[int(kmag)] = Efk[int(kmag)] + 0.5 * fk[ii] * np.conj(fk[ii])

   figE, axE = plt.subplots(1,1, figsize=(10,0.75*10))
   axE.plot(k, np.real(Efk), label=r'$E_{f}$')
   axE.set_xscale('log')
   axE.set_yscale('log')
   #axE.set_xlim(k[1], k.max())
   #axE.set_ylim(1.0e-17, max(np.real(ETk).max(), np.real(EKk).max()))
   #axE.set_title(r'Energy spectrum of $\frac{1}{2}\sin^{2}\left(2x\right)$')
   axE.set_xlabel(r'$k$')
   axE.set_ylabel(r'$E\left(k\right)$')
   axE.legend()

   figE.tight_layout()
   figE.savefig('simple_test.pdf')

   plt.show()


if __name__ == '__main__':
    plot_energy_spectra()
