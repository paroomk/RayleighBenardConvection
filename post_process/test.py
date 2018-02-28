#!/usr/bin/python

import sys
import h5py
import math
import config    
import numpy as np 
import numpy.matlib 
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import matplotlib.cm as cm

def plot_energy_spectra():   # Note: computes thermal energy from perturbation temperature

   #file = path + filename + str(8) + '.h5'
   #f = h5py.File(file,'r+')

   #x = f.get('scales/x/1.0')
   #x = np.array(x)
   
   #z = f.get('scales/z/1.0')
   #z = np.array(z)

   #Nx = x.shape[0]

   #z1 = -0.8
   #jj = (np.abs(z-z1)).argmin()
   
   #T = f.get('tasks/b')
   #T  = np.array(T)

   #T = T[1,:,jj]

   dN = 1 
   
   Nx = 64
  
   x = np.arange(Nx)/Nx * (2*np.pi) 
   
   T_k  = np.zeros([Nx], dtype=np.complex_)
   
   T = np.cos(3*x)*np.sin(x) + 4*np.cos(x) + 5

   fig1, ax1 = plt.subplots(1,figsize=[7,4])
   ax1.plot(x,T, 'g', label = 'Input function')
   ax1.set_xlabel('x')
   ax1.set_ylabel('T')
   ax1.set_title('test case')

   plt.show()

   T_k = np.fft.fft(T)/(Nx)
           
   ThE_bk  =  0.5*T_k*np.conj(T_k)/dN
        
   print(ThE_bk)
   
   k = np.concatenate((np.arange(0,Nx//2+1), np.arange(-Nx//2 + 1,0)))  

   ThE_bk = [x for _,x in sorted(zip(k,ThE_bk))]

   ThE_k  = np.zeros([Nx//2+1], dtype=np.complex_)

   ThE_k[1:Nx//2-1] = np.add(ThE_bk[Nx//2-2:0:-1],ThE_bk[Nx//2:Nx-2])

   ThE_k[Nx//2] = 2*ThE_bk[Nx-1]

   ThE_k[0] = ThE_bk[Nx//2-1]

   k = np.arange(0,Nx//2+1)

   alpha = 1
   
   freq = k*alpha

   print(ThE_k)
   
   fig2, ax2 = plt.subplots(1,figsize=[7,4])
   ax2.plot(freq,np.real(ThE_k), 'b*', label = 'Thermal Energy')
   ax2.set_xlabel('k')
   ax2.set_ylabel('E_k')
   ax2.set_title('test case')
#   ax2.set_ylim((1E-9,1E-1))
   ax2.legend(loc='upper right')
  
   plt.tight_layout()
   plt.show()

   return 

path = '/Volumes/Chest/snapshots_1E7_10_24_2/'
filename = 'snapshots_1E7_10_24_2_s'
plot_energy_spectra()

