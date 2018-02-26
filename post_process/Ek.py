#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt
import h5py
import config    

plt.rcParams.update(config.pars)

def plot_energy_spectra():

   # Read file 
   fname = 'snapshots_1E7_4_32_2_s20.h5'
   f = h5py.File(fname, 'r')

   # Coordinates
   x = np.array(f['scales/x/1.0'])
   z = np.array(f['scales/z/1.0'])

   Nx = x.shape[0]
   Nz = z.shape[0]

   # Solution fields
   T = np.array(f['tasks/b'])
   #T = T[10,:,Nz//2] # Get T near center and at specific time
   T = T[1,:,50] # Get T near the wall and at specific time
   #Ta = T - z[50]

   # Energy
   ET = 0.5 * T * T

   # Plot in physical space
   fig1, ax1 = plt.subplots(1,1)
   ax1.plot(x,T)

   fig2, ax2 = plt.subplots(1,1, figsize=(10,0.75*10))
   ax2.plot(x,ET)
   ax2.set_ylabel(r'$K_{T} = \frac{1}{2}T^{2}$')
   ax2.set_xlabel(r'$x$')
   ax2.set_xlabel(r'$x$')
   fig2.tight_layout()
   fig2.savefig('K_T_ymid.pdf')

   # Take FFT
   Tk = np.fft.fft(T) / Nx
   ETk = np.fft.fft(ET) / Nx

   # Set up wavenumbers
   kx = np.zeros(Nx)
   for ii in range(Nx//2):
      kx[ii] = float(ii)
   for ii in range(Nx//2, Nx):
      kx[ii] = float(ii-Nx)

   # Compute energy in Fourier space
   NF = Nx
   kmax = NF // 2
   Ek = np.zeros(NF//2, dtype=np.complex_)
   k  = np.arange(0.0, kmax)
   for ii in range(Nx):
      kmag = abs(kx[ii])
      if (kmag < kmax):
         Ek[int(kmag)] = Ek[int(kmag)] + 0.5*Tk[ii] * np.conj(Tk[ii])

   figE, axE = plt.subplots(1,1, figsize=(10,0.75*10))
   axE.plot(k, np.real(Ek))
   axE.set_xscale('log')
   axE.set_yscale('log')
   axE.set_xlabel(r'$k$')
   axE.set_ylabel(r'$E_{T}\left(k\right)$')
   axE.set_title(r'Spectrum at $y\approx y_{\textrm{wall}}$')
   #axE.set_title(r'Spectrum at $y\approx y_{\textrm{mid}}$')

   figE.tight_layout()
   figE.savefig('ET_ywall.pdf')

   # Try to average over neighboring shells
   Ek_ave = np.zeros(NF//2)
   for idx, val in enumerate(k):
       if val == 0:
          Ek_ave[idx] = np.real(0.5 * (Ek[idx] + Ek[idx+1]))
       elif val == kmax - 1:
          Ek_ave[idx] = np.real(0.5 * (Ek[idx-1] + Ek[idx]))
       else:
          Ek_ave[idx] = np.real((Ek[idx-1] + Ek[idx] + Ek[idx+1]) / 3.0)

   figA, axA = plt.subplots(1,1, figsize=(10,0.75*10))
   axA.plot(k, Ek_ave)
   axA.set_xscale('log')
   axA.set_yscale('log')
   axA.set_xlabel(r'$k$')
   axA.set_ylabel(r'$\left<E_{T}\left(k\right)\right>$')
   axA.set_title(r'Shell-averaged Spectrum at $y\approx y_{\textrm{wall}}$')

   figA.tight_layout()
   figA.savefig('ET_ave_ywall.pdf')

   plt.show()


if __name__ == '__main__':
    plot_energy_spectra()