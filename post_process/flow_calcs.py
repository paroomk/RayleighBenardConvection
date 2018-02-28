#!/usr/bin/python

import numpy as np 
import h5py
import config    

np.set_printoptions(threshold=np.inf)

def calc_Re():

   # Set up some parameters
   Ra = 1.0e+07
   Pr = 4.0

   nu = np.sqrt(16.0 * Pr / Ra)

   # Read file 
   fname = 'snapshots_1E7_4_32_2_s20.h5'
   f = h5py.File(fname, 'r')

   # Coordinates
   x = np.array(f['scales/x/1.0'])
   z = np.array(f['scales/z/1.0'])
   dz = np.ediff1d(z)
   Lz = 2.0
   h = Lz / 2.0

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

   #Ta = T - z[50]

   # Compute RMS of velocity field
   ubar = 0.0
   u2bar = 0.0
   vbar = 0.0
   v2bar = 0.0
   for ii in range(Nx):
       for jdx, deltaz in enumerate(dz):
           ubar += uxy[ii,jdx+1] * deltaz
           u2bar += uxy[ii,jdx+1] * uxy[ii,jdx+1] * deltaz

           vbar += vxy[ii,jdx+1] * deltaz
           v2bar += vxy[ii,jdx+1] * vxy[ii,jdx+1] * deltaz

   ubar = ubar / Lz / Nx
   u2bar = u2bar / Lz / Nx

   urms2 = u2bar - ubar * ubar
   vrms2 = v2bar - vbar * vbar

   urms = np.sqrt(urms2)
   vrms = np.sqrt(vrms2)
   rms = np.sqrt(0.5 * (urms2 + vrms2))

   # Compute Re based on different rms velocities
   Re = rms * h / nu
   Re_u = urms * h / nu
   Re_v = vrms * h / nu

   print("Re_rms = {0:25.16e},   rms = {1:25.16e}".format(Re, rms))
   print("Re_urms = {0:25.16e},   urms = {1:25.16e}".format(Re_u, urms))
   print("Re_vrms = {0:25.16e},   vrms = {1:25.16e}".format(Re_v, vrms))

   print()
   # Compute Re based off of max velocity in domain
   umax = abs(uxy).max()
   vmax = abs(vxy).max()
   velmax = max(umax, vmax)

   Re_umax = umax * h / nu
   Re_vmax = vmax * h / nu
   Re_vmax = velmax * h / nu

   print("Re_umax = {0:25.16e},   umax = {1:25.16e}".format(Re_umax, umax))
   print("Re_vmax = {0:25.16e},   vmax = {1:25.16e}".format(Re_umax, vmax))
   print("Re_max = {0:25.16e},   velmax = {1:25.16e}".format(Re_umax, velmax))

if __name__ == '__main__':
    calc_Re()
