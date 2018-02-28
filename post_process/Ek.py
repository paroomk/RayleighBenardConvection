#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import ticker
import h5py
import config    

plt.rcParams.update(config.pars)
np.set_printoptions(threshold=np.inf)

def plot_energy_spectra():

   # Read file 
   fname = 'snapshots_1E7_4_32_2_s20.h5'
   f = h5py.File(fname, 'r')

   # Coordinates
   x = np.array(f['scales/x/1.0'])
   z = np.array(f['scales/z/1.0'])
   dz = np.ediff1d(z)
   Lz = 2.0

   Nx = x.shape[0]
   Nz = z.shape[0]

   # Solution fields
   T = np.array(f['tasks/b'])
   u = np.array(f['tasks/u'])
   v = np.array(f['tasks/w'])

   # Get fields at specific time
   Txy = T[1,:,:]
   uxy = u[1,:,:]
   vxy = v[1,:,:]
   #T = T[10,:,Nz//2] # Get T near center and at specific time
   T = T[1,:,50] # Get T near the wall and at specific time
   u = u[1,:,50] # Get u near the wall and at specific time
   v = v[1,:,50] # Get v near the wall and at specific time

   #Ta = T - z[50]

   # Take FFT
   Tk = np.fft.fft(T) / Nx
   uk = np.fft.fft(u) / Nx
   vk = np.fft.fft(v) / Nx

   # Compute FFT of fields in x
   Tky = np.fft.fft(Txy, axis=0) / Nx
   uky = np.fft.fft(uxy, axis=0) / Nx
   vky = np.fft.fft(vxy, axis=0) / Nx
   
   # Compute y-average
   #Tbar = np.zeros(Nx)
   #ubar = np.zeros(Nx)
   #vbar = np.zeros(Nx)
   #for ii in range(Nx):
   #    for idy, deltaz in enumerate(dz):
   #        Tbar[ii] += Txy[ii,idy+1] * deltaz
   #        ubar[ii] += uxy[ii,idy+1] * deltaz
   #        vbar[ii] += vxy[ii,idy+1] * deltaz
   #Tbar /= Lz
   #ubar /= Lz
   #vbar /= Lz

   # Take FFT of y-averages
   #Tbark = np.fft.fft(Tbar) / Nx
   #ubark = np.fft.fft(ubar) / Nx
   #vbark = np.fft.fft(vbar) / Nx

   # Set up wavenumbers
   kx = np.zeros(Nx)
   for ii in range(Nx//2):
      kx[ii] = float(ii)
   for ii in range(Nx//2, Nx):
      kx[ii] = float(ii-Nx)

   # Compute energy in Fourier space
   NF = Nx
   kmax = NF // 2
   #kmax = np.ceil(NF / 2)
   #NF = int(np.floor(2.0 * Nx / 3.0))
   ETk = np.zeros(int(kmax), dtype=np.complex_)
   ETbark = np.zeros(int(kmax), dtype=np.complex_)
   EKk = np.zeros(int(kmax), dtype=np.complex_)
   EKbark = np.zeros(int(kmax), dtype=np.complex_)
   k  = np.arange(0.0, kmax)
   for ii in range(Nx):
      kmag = abs(kx[ii])

      Tbar = 0.0
      ubar = 0.0
      vbar = 0.0
      for idy, deltaz in enumerate(dz):
          Tbar += np.real(Tky[ii, idy+1] * np.conj(Tky[ii, idy+1])) * deltaz
          ubar += np.real(uky[ii, idy+1] * np.conj(uky[ii, idy+1])) * deltaz
          vbar += np.real(vky[ii, idy+1] * np.conj(vky[ii, idy+1])) * deltaz
      Tbar = Tbar / Lz
      ubar = ubar / Lz
      vbar = vbar / Lz

      if (kmag < kmax):
         ETk[int(kmag)] = ETk[int(kmag)] + 0.5 * Tk[ii] * np.conj(Tk[ii])
         ETbark[int(kmag)] = ETbark[int(kmag)] + 0.5 * Tbar

         EKk[int(kmag)] = EKk[int(kmag)] + \
             0.5 * (uk[ii] * np.conj(uk[ii]) + vk[ii] * np.conj(vk[ii]))
         EKbark[int(kmag)] = EKbark[int(kmag)] + 0.5 * (ubar + vbar)

   figE, axE = plt.subplots(1,1, figsize=(10,0.75*10))
   axE.plot(k, np.real(ETk), label=r'$E_{T}$')
   axE.plot(k, np.real(EKk), label=r'$E_{K}$')
   axE.set_xscale('log')
   axE.set_yscale('log')
   axE.set_xlim(k[1], k.max())
   axE.set_ylim(1.0e-17, max(np.real(ETk).max(), np.real(EKk).max()))
   axE.set_xlabel(r'$k$')
   axE.set_ylabel(r'$E\left(k\right)$')
   axE.set_title(r'Spectrum at $y\approx y_{\textrm{wall}}$')
   #axE.set_title(r'Spectrum at $y\approx y_{\textrm{mid}}$')
   axE.legend()

   figE.tight_layout()
   figE.savefig('Ek_ywall.pdf')

   figEbar, axEbar = plt.subplots(1,1, figsize=(10,0.75*10))
   axEbar.plot(k, np.real(ETbark), label=r'$\overline{E}_{T}$')
   axEbar.plot(k, np.real(EKbark), label=r'$\overline{E}_{K}$')
   axEbar.set_xscale('log')
   axEbar.set_yscale('log')
   axEbar.set_xlim(k[1], k.max())
   axEbar.set_ylim(1.0e-17, max(np.real(ETbark).max(), np.real(EKbark).max()))
   axEbar.set_xlabel(r'$k$')
   axEbar.set_ylabel(r'$\overline{E}\left(k\right)$')
   axEbar.set_title(r'Spectrum at $y\approx y_{\textrm{wall}}$')
   #axE.set_title(r'Spectrum at $y\approx y_{\textrm{mid}}$')
   axEbar.legend()

   figEbar.tight_layout()
   figEbar.savefig('Ekbar_ywall.pdf')


   # Try to average over neighboring shells
   ETk_ave = np.zeros(int(kmax))
   EKk_ave = np.zeros(int(kmax))
   for idx, val in enumerate(k):
       if val == 0:
          ETk_ave[idx] = np.real(0.5 * (ETk[idx] + ETk[idx+1]))
          EKk_ave[idx] = np.real(0.5 * (EKk[idx] + EKk[idx+1]))
       elif val == kmax - 1:
          ETk_ave[idx] = np.real(0.5 * (ETk[idx-1] + ETk[idx]))
          EKk_ave[idx] = np.real(0.5 * (EKk[idx-1] + EKk[idx]))
       else:
          ETk_ave[idx] = np.real((ETk[idx-1] + ETk[idx] + ETk[idx+1]) / 3.0)
          EKk_ave[idx] = np.real((EKk[idx-1] + EKk[idx] + EKk[idx+1]) / 3.0)

   figA, axA = plt.subplots(1,1, figsize=(10,0.75*10))
   axA.plot(k, ETk_ave, label=r'$\left<E_{T}\left(k\right)\right>$')
   axA.plot(k, EKk_ave, label=r'$\left<E_{K}\left(k\right)\right>$')
   axA.set_xscale('log')
   axA.set_yscale('log')
   axA.set_xlim(k[1], k.max())
   axA.set_ylim(1.0e-17, max(ETk_ave.max(), EKk_ave.max()))
   axA.set_xlabel(r'$k$')
   axA.set_ylabel(r'$\left<E\left(k\right)\right>$')
   axA.set_title(r'Shell-averaged Spectrum at $y\approx y_{\textrm{wall}}$')
   axA.legend()

   figA.tight_layout()
   figA.savefig('Ek_shellave_ywall.pdf')

   # Create energy contours in Fourier space
   NF = Nx
   #kmax = NF // 2
   #Eky = np.zeros([NF//2, Nz], dtype=np.complex_)
   ETky = np.zeros([int(kmax), Nz], dtype=np.complex_)
   K, Z = np.meshgrid(k, z)
   for jj in range(Nz):
       for ii in range(Nx):
           kmag = abs(kx[ii])
           if (kmag < kmax):
              ETky[int(kmag), jj] = ETky[int(kmag), jj] + \
                  0.5 * Tky[ii,jj] * np.conj(Tky[ii,jj])

   figC, axC = plt.subplots(1,1, figsize=(10,0.75*10))
   CS = axC.contourf(K, Z, np.real(np.transpose(ETky)), locator=ticker.LogLocator())
   figC.colorbar(CS);

   # Create contours in physical space
   X, Z = np.meshgrid(x, z)

   figCp, axCp = plt.subplots(1,1, figsize=(10,0.75*10))
   CSp = axCp.contourf(X, Z, np.transpose(Txy), 101, cmap='jet')
   figCp.colorbar(CSp);

   plt.show()


if __name__ == '__main__':
    plot_energy_spectra()
