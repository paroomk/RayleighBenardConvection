#!/usr/bin/python

import sys
import h5py
import math
import config    
import numpy as np
from numpy import linalg as LA
from scipy import interpolate 
from scipy import integrate 
from detect_peaks import detect_peaks 
import numpy.matlib 
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import matplotlib.cm as cm

def integrate_z(Fz,z):

    F = np.trapz(Fz, z, axis=0)
    
    return F

def grad(F,x,z):
    dx = x[1]-x[0]

    Nx = x.shape[0]
    Ny = z.shape[0]

    Fx = np.zeros((Nx,Ny))
    Fy = np.zeros((Nx,Ny))

    for i,xi in enumerate(x):
        Fx[i,:] = (F[(i+1) % Nx,:] - F[i-1,:])/(2*dx)

    for i,zi in enumerate(z):
        if i == 0:

            dx1 = z[i+1] - z[i] 
            dx2 = z[i+2] - z[i+1]

            a = -(2*dx1 + dx2)/(dx1*(dx1+dx2))
            b = (dx1 + dx2)/(dx1*dx2)
            c = - dx1/(dx2*(dx1+dx2))

            Fy[:,i] = a*F[:,i] + b*F[:,i+1] + c*F[:,i+2]

        elif i == Ny-1:

            dx1 = z[i] - z[i-1]              # N
            dx2 = z[i-1] - z[i-2]           # N-1

            a = (2*dx1 + dx2)/(dx1*(dx1+dx2))
            b = -(dx1 + dx2)/(dx1*dx2)
            c = dx1/(dx2*(dx1+dx2))

            Fy[:,i] = a*F[:,i] + b*F[:,i-1] + c*F[:,i-2]

        else :

            dx1 = z[i] - z[i-1]
            dx2 = z[i+1] - z[i] 

            a = -dx2/(dx1*(dx1+dx2))
            b = (dx2-dx1)/(dx1*dx2)
            c = dx1/(dx2*(dx1+dx2))

            Fy[:,i] = a*F[:,i-1] + b*F[:,i] + c*F[:,i+1]
    
    #print(-np.sum(Fy[:,0],axis = 0)/Nx )

    return Fx, Fy

def opt_comparison(): #Computes L2 of error between turbulent plumes and optimal structure

    file = path + filename + str(520) + '.h5'
    f = h5py.File(file,'r+')
   
    x = f.get('scales/x/1.0')
    x = np.array(x)
   
    z = f.get('scales/z/1.0')
    z = np.array(z)

    T = f.get('tasks/b')
    T = np.array(T)
   
    ti = f.get('scales/sim_time')
    t = np.array(ti)
   
    Nx = T.shape[1]
    Ny = T.shape[2]

    dx = x[1] - x[0]

    for i in range(520, 521):
        file = path + filename + str(i) + '.h5'
        print(file)
        f = h5py.File(file,'r+')
   
        Ti = f.get('tasks/b')
        Ti = np.array(Ti)

        ti = f.get('scales/sim_time')
    
        t  = np.concatenate((t, np.array(ti)),0)
   
        T  = np.concatenate((T, np.array(Ti)),0)

    Mx = 128 ; My = 201;

    Lx = 2*np.pi/alpha;
    x_o = np.linspace(-Lx/2, Lx/2- Lx/Mx, Mx)
     
    y_l = Ny//2  # optimal solution is aligned based on the plume behavior along z[y_l] 
    loc = Lx/2

    i_r = np.nonzero((x > (-loc - dx/2)) & (x < (loc + dx/2)))  # range of indices covered by the optimal solution when placed at the x = 0
    N0  = np.nonzero(x == 0.)                                    #index of the center of the box

    T_o = np.loadtxt('/Volumes/Work/Fluids_project/Programs/POD/optimal_solns/optimal_solutions/T_1E7_10.txt')

    c =  (np.arange(0,My)) / (My - 1.)
    y = np.cos(c*np.pi);

    T_o = np.reshape(T_o,[Mx,My])
    [T_ox,T_oy] = grad(T_o,x_o,-y)


    x_i = x[i_r]

    Xi,Yi = np.meshgrid(x_i,z)

    fT_i  = interpolate.interp2d(x_o,-y,T_o.T,kind = 'quintic')
    fT_ix = interpolate.interp2d(x_o,-y,T_ox.T,kind = 'quintic')
    fT_iy = interpolate.interp2d(x_o,-y,T_oy.T,kind = 'quintic')

    T_i  = fT_i(x_i,z) 
    T_ix = fT_ix(x_i,z) 
    T_iy = fT_iy(x_i,z)
 
    E_l2  = np.zeros(len(t))  #L2 norm of error
    E_h1  = np.zeros(len(t))  #H1 semi-norm of error
    x_min  = np.zeros(len(t))  #H1 semi-norm of error
    E_l2_j= np.zeros(Nx)              
    E_h1_j= np.zeros(Nx)              #10 = max expected #plume detected
    #group = np.zeros(10)
    #group[:] = 100

    for i in range(0,len(t)):
        T_s = T[i,:,:] - np.matlib.repmat(z,Nx,1)
        [T_sx,T_sy] = grad(T_s,x,z)
        T_y = T_s[:,y_l] 
#        print(T_sx,T_sy)

        X,Y = np.meshgrid(x_o,-y)
#        X,Y = np.meshgrid(x,z)

 
#        Check interpolated plume here
           
        fig, axarr = plt.subplots(1, figsize=(6,7))
        levels = np.linspace(-1,1,11)
        axarr.contourf(X, Y, T_o.T,levels,cmap=cm.seismic)
#        axarr.contourf(X,Y, T_sy.T,cmap=cm.seismic)
        plt.show()
        quit()

#        loc = x[plume_locs]

        E_l2_j[:] = 0
        E_h1_j[:] = 0

#        print(len(loc))

        for j in range(0,Nx):

            shift = np.subtract(j,N0) #index of plume location from center
            j_r = np.add(i_r,shift)              

            T_p = T_s[j_r % Nx,:].T 
 
            T_px = T_sx[j_r % Nx,:].T  
            T_py = T_sy[j_r % Nx,:].T  

            #E_l2_j[j] = LA.norm(T_i-T_p[:,:,0])

            E_l2_j[j] = np.sqrt(integrate_z(np.mean((T_i-T_p[:,:,0])**2,axis = 1),z)*Lx)
            #E_h1_j[j] = LA.norm(np.dstack((T_ix-T_px[:,:,0],T_iy-T_py[:,:,0]))) 
            E_h1_j[j] = np.sqrt(integrate_z(np.mean((T_ix-T_px[:,:,0])**2,axis = 1),z)*Lx + integrate_z(np.sum((T_iy-T_py[:,:,0])**2,axis = 1),z)*Lx)
            #print(E_l2_j[j])

        E_l2[i] = min(E_l2_j)
        j_min = np.argmin(E_l2_j)
        x_min[i] = x[j_min]

        E_h1[i] = min(E_h1_j)


    #print(E_l2)
    fig, axarr = plt.subplots(2, figsize=(6,7))
    axarr[0].plot(t[E_l2[:]>0],E_l2[E_l2[:]>0])
    axarr[0].legend(loc='upper right')
    axarr[0].set_xlabel('time')
    axarr[0].set_ylabel('E_l2')

    axarr[1].plot(t[E_h1[:]>0],E_h1[E_h1[:]>0])
    #axarr[1].plot(t,x_min,'*')
    axarr[1].legend(loc='upper right')
    axarr[1].set_xlabel('time')
    axarr[1].set_ylabel('E_h1')
    plt.show()

    return
    
def plot_contours():   #plots the last snapshot generated
    file = path + filename + str(nfiles) + '.h5'
    print(file)
    f = h5py.File(file,'r+')
        
    x = f.get('scales/x/1.0')
    x = np.array(x)
   
    z = f.get('scales/z/1.0')
    z = np.array(z)

    t = f.get('scales/sim_time')
    t = np.array(t)

    N = t.shape[0]
    Nx = x.shape[0]

    X,Z = np.meshgrid(x,z)

    T = f.get('tasks/b')
    T = np.array(T)
    
    u = f.get('tasks/u')
    u = np.array(u)
    
    w = f.get('tasks/w')
    w = np.array(w)
    
    fig, axarr = plt.subplots(3, figsize=(6,7))
    levels = np.linspace(-1,1,11)
    axarr[0].contourf(X, Z, T[N-1,:,:].T - np.matlib.repmat(z,Nx,1).T,levels,cmap=cm.seismic)
    axarr[0].set_title('Temperature at t = %s'%(t[N-1]))
    axarr[1].contourf(X, Z, u[N-1,:,:].T,cmap=cm.seismic)
    axarr[1].set_title('Horizontal Velocity at t = %s'%(t[N-1]))
    axarr[2].contourf(X, Z, w[N-1,:,:].T,cmap=cm.seismic)
    axarr[2].set_title('Vertical Velocity at t = %s'%(t[N-1]))
    plt.tight_layout()
    plt.show()
    return

def plot_energyspectra():
    f = open('fileread', 'r')
    #fig1, ax1 = plt.subplots(1,figsize = [5,4])
    fig, axarr = plt.subplots(2, figsize=(6,7))
    kmax = [] 
    for i, line in enumerate(f):
        print(line)
        line = line.strip('\n')
        E = np.loadtxt(E_type + line + '.txt')
        kx = np.loadtxt('kx' + line + '.txt')
        axarr[0].loglog(kx,E, label = line)
        index = E.argmax()
        kmax.append(kx[index])
    axarr[0].set_xlabel('kx')
    axarr[0].set_ylabel(E_type)
    axarr[0].legend(bbox_to_anchor=(0.5, 0.5), loc='lower left')
    axarr[0].set_ylim((1E-15,1))

    axarr[1].plot(size,kmax)
    axarr[1].set_xlabel('box size')
    axarr[1].set_ylabel('kx_max')
    plt.show()
    return
 
def plot_Nus(): 
    f = open('fileread', 'r')
    fig1, ax1 = plt.subplots(1,figsize = [5,4])
    for line in f:
        print(line)
        line = line.strip('\n')
        Nu = np.loadtxt('Nu' + line + '.txt')
        t = np.loadtxt('t' + line + '.txt')
        ax1.plot(t[200:1050],Nu[200:1050], label = line)
    ax1.set_xlabel('time')
    ax1.set_ylabel('Nu')
    ax1.legend(bbox_to_anchor=(0.5, 0.5), loc='lower right')
    plt.show()
    
    return 
  
def plot_Nu(): 
   file = path + filename + str(50) + '.h5'
   f = h5py.File(file,'r+')
   
   dT = f.get('tasks/bz')
   dT = np.array(dT)
   
   ti = f.get('scales/sim_time')
   t = np.array(ti)
   
   Nx = dT.shape[1]
   Ny = dT.shape[2]
   print(Nx,Ny)

   Nu = - dT[:,:,1].sum(1)/Nx + 1; 

   for x in range(51, 79):
       file = path + filename + str(x) + '.h5'
       print(file)
       f = h5py.File(file,'r+')
   
       dT = f.get('tasks/bz')
       dT = np.array(dT)

       ti = f.get('scales/sim_time')

       Nui = - dT[:,:,1].sum(1)/Nx + 1; 
   
       t  = np.concatenate((t, np.array(ti)),0)
   
       Nu  = np.concatenate((Nu, np.array(Nui)),0)

   t1 = 200 
   t2 = 300
  
   i1 = (np.abs(t-t1)).argmin()
   i2 = (np.abs(t-t2)).argmin()

   Nu_avg = Nu[i1:i2].sum(0)/(i2-i1+1)   
   
   print('Time-averaged Nusselt number = %s '%(Nu_avg))

   fig1, ax1 = plt.subplots(1,figsize = [5,4])
   ax1.plot(t,Nu)
   ax1.set_xlabel('time')
   ax1.set_ylabel('Nu')
   ax1.set_title(case)

   print(sys.getsizeof(Nu_avg))
   print(Nu.shape[0])

   np.savetxt("Nu"+ caseid + ".txt", Nu)
   np.savetxt("t"+ caseid + ".txt", t)

   plt.show()
   return

def track_localNu(xloc): #finds Nu(t) in an optimal box centered at xloc
    
   file = path + filename + str(50) + '.h5'
   f = h5py.File(file,'r+')
   
   dT = f.get('tasks/bz')
   dT = np.array(dT)
   
   ti = f.get('scales/sim_time')
   t = np.array(ti)
   
   x = f.get('scales/x/1.0')
   x = np.array(x)
   
   Nx = dT.shape[1]
   Ny = dT.shape[2]
   
   Lx = 2*np.pi/alpha;
   
   dx = x[1]-x[0]
   
   loc = Lx/8
   
   x = x - dx*Nx/2
   
   i_r = np.nonzero((x > (-loc - dx/2)) & (x < (loc + dx/2)))  # range of indices covered by the optimal solution when placed at the x = 0
   
   N0 = np.argmin(np.abs(x-0.0))
   j  = np.argmin(np.abs(x-xloc))
   
   shift = np.subtract(j,N0) #index of plume location from center
   i_opt = np.add(i_r,shift)
   
   print(np.argmin(np.abs(x-xloc)))
   
   print(x[i_opt % Nx])
   
   Nx_opt = i_opt.shape[1]
   
   print(Nx_opt)

   Nu = - dT[:,i_opt % Nx,Ny-1].sum(2)/Nx_opt + 1;
   
   #print(np.shape(dT[:,i_opt % Nx,1]))
   
   #print(Nu) 

   for i in range(51, 79):
       file = path + filename + str(i) + '.h5'
       print(file)
       f = h5py.File(file,'r+')
   
       dT = f.get('tasks/bz')
       dT = np.array(dT)

       ti = f.get('scales/sim_time')

       Nui = - dT[:,i_opt % Nx,Ny-1].sum(2)/Nx_opt + 1; 
   
       t  = np.concatenate((t, np.array(ti)),0)
   
       Nu  = np.concatenate((Nu, np.array(Nui)),0)

   t1 = 600 
   t2 = 900
  
   i1 = (np.abs(t-t1)).argmin()
   i2 = (np.abs(t-t2)).argmin()

   Nu_avg = Nu[i1:i2].sum(0)/(i2-i1+1)   
   
   #print('Time-averaged Nusselt number = %s '%(Nu_avg))

   fig1, ax1 = plt.subplots(1,figsize = [5,4])
   ax1.plot(t,Nu)
   ax1.set_xlabel('time')
   ax1.set_ylabel('Nu')
   ax1.set_title(case)

   print(sys.getsizeof(Nu_avg))
   print(Nu.shape[0])

   plt.show()
    
    
   return

def spatial_correlation():
    file = path + filename + str(50) + '.h5'
    f = h5py.File(file,'r+')
   
    x = f.get('scales/x/1.0')
    x = np.array(x)
   
    z = f.get('scales/z/1.0')
    z = np.array(z) - 1

    T = f.get('tasks/b')
    T = np.array(T)
   
    ti = f.get('scales/sim_time')
    t = np.array(ti)
   
    Nx = T.shape[1]
    Ny = T.shape[2]

    dx = x[1] - x[0]
    
    x = x - dx*Nx/2
    
    #print(x)

    for i in range(51,79):
        file = path + filename + str(i) + '.h5'
        print(file)
        f = h5py.File(file,'r+')
   
        Ti = f.get('tasks/b')
        Ti = np.array(Ti)

        ti = f.get('scales/sim_time')
    
        t  = np.concatenate((t, np.array(ti)),0)
   
        T  = np.concatenate((T, np.array(Ti)),0)

    Mx = 128 ; My = 201;

    Lx = 2*np.pi/alpha;
    x_o = np.linspace(-Lx/2, Lx/2 - Lx/Mx, Mx)
    
    print(x_o)
     
    #y_l = Ny//2  # optimal solution is aligned based on the plume behavior along z[y_l] 
    loc = Lx/8 #adjust loc to restrict region of interest.

    i_r = np.nonzero((x > (-loc - dx/2)) & (x < (loc + dx/2)))   #range of indices covered by the optimal solution when placed at the x = 0
    print(x[i_r])
    N0  = np.nonzero(x == 0.)                                    #index of the center of the box

    T_o = np.loadtxt('/Volumes/Work/Fluids_project/Programs/POD/optimal_solns/optimal_solutions/T_5E6_7.txt')

    c =  (np.arange(0,My)) / (My - 1.)
    y = np.cos(c*np.pi);
    
    #print(z)

    T_o = np.reshape(T_o,[Mx,My])


    x_i = x[i_r]

    Xi,Yi = np.meshgrid(x_i,z)

    fT_i  = interpolate.interp2d(x_o,-y,T_o.T,kind = 'linear')

    T_i  = fT_i(x_i,z) 
    
    nT_i =  np.sqrt(integrate_z(np.mean(T_i**2,axis = 1)*Lx,z))  #int (T_i dx dz) in the optimal box  
 
    cT  = np.zeros(len(t))  #L2 norm of error
    x_max  = np.zeros(len(t))  #H1 semi-norm of error
    cT_j= np.zeros(Nx) 
    
    X,Y = np.meshgrid(x_i,z)
    
    #Check interpolated plume here
    
    fig, axarr = plt.subplots(1, figsize=(6,7))
    levels = np.linspace(-1,1,11)
    axarr.contourf(X, Y, T_i,levels,cmap=cm.seismic)
    plt.show()             

    quit()

    for i in range(0,len(t)):
        T_s = T[i,:,:] - np.matlib.repmat(z,Nx,1)
        #T_y = T_s[:,y_l] 

        cT_j[:] = 0

        for j in range(0,Nx):

            shift = np.subtract(j,N0) #index of plume location from center
            j_r = np.add(i_r,shift)              

            T_p = T_s[j_r % Nx,:].T 
 
            #cT_j[j] = LA.norm(T_i-T_p[:,:,0])
            
            nT_p =  np.sqrt(integrate_z(np.mean(T_p[:,:,0]**2,axis = 1)*Lx,z))    #int (T_p dx dz) in the optimal box                      

            cT_j[j] = (integrate_z(np.mean((T_i*T_p[:,:,0]),axis = 1)*Lx,z))/(nT_i*nT_p)
            #print(cT_j[j])

        cT[i] = max(cT_j)
        j_max = np.argmax(cT_j)
        x_max[i] = x[j_max]

    #print(cT)
    fig, axarr = plt.subplots(2, figsize=(6,7))
    axarr[0].plot(t[cT[:]>0],cT[cT[:]>0])
    axarr[0].legend(loc='upper right')
    axarr[0].set_xlabel('time')
    axarr[0].set_ylabel('Correlation of temperature fields')
    
    axarr[1].plot(t[cT[:]>0],x_max[cT[:]>0])
    axarr[1].legend(loc='upper right')
    axarr[1].set_xlabel('time')
    axarr[1].set_ylabel('x_max')
    plt.show()
    
    cT_max = max(cT)
    j_max  = np.argmax(cT)
    x_gmax = x_max[j_max]
    
    print(cT_max,t[j_max],x_gmax)

    return

def plot_BLs(): 
    f = open('fileread', 'r')
    fig1, ax1 = plt.subplots(1,figsize = [5,4])
    for line in f:
        print(line)
        line = line.strip('\n')
        b = np.loadtxt('BL' + line + '.txt')
        z_H = np.loadtxt('z_H' + line + '.txt')
        ax1.semilogx(z_H,b, label = line)
    ax1.set_xlabel('z/H')
    ax1.set_ylabel('(u_x)^2 + (w_z)^2')
    ax1.legend(bbox_to_anchor=(0.5, 0.5), loc='lower right')
    plt.show()
    
    return 

def plot_optTBL():    #plots momentum boundary layer through x[j]
   T_o = np.loadtxt('/Volumes/Work/Fluids_project/Programs/POD/optimal_solns/optimal_solutions/T_5E6_7.txt')

   Mx = 128
   My = 201

   Lx = 2*np.pi/alpha;

   c =  (np.arange(0,My)) / (My - 1.)
   y = -np.cos(c*np.pi) + 1
   x_o = np.linspace(-Lx/2, Lx/2- Lx/Mx, Mx)

   T_o = np.reshape(T_o,[Mx,My])

   X,Y = np.meshgrid(x_o,y)
           
   #fig, axarr = plt.subplots(1, figsize=(6,7))
   #p1 = axarr.contourf(X, Y, T_o.T,cmap=cm.seismic)
   #cbar = plt.colorbar(p1)
   #plt.show()
   
   b = 0
   j = 0 #Mx//2 #Only Centre

   for i,xi in enumerate(x_o[j:j+1]):
       Ti = T_o[i,:]
       b =  b + Ti#/Mx

   H = 2

   y_s = 1E-3*H
   y_e = 5E-1*H

   j_s = (np.abs(y-y_s)).argmin()
   j_e = (np.abs(y-y_e)).argmin()
       
   j_max = b[j_s:j_e].argmax()
#   print(y[j_max]/H)

   np.savetxt("BL_T"+ optcaseid + ".txt", b[j_s:j_e])
   np.savetxt("z_H"+ optcaseid + ".txt", y[j_s:j_e]/H)

   fig1, ax1 = plt.subplots(1,figsize = [5,4])
   ax1.semilogx(y[j_s:j_e]/H,b[j_s:j_e])
   ax1.set_xlabel('z/H')
   ax1.set_ylabel('T(x_0,z)')
   ax1.set_title(case)
   plt.show()

   return

def find_optMBL():  #finds locations and times of optimal match
   file = path + filename + str(79) + '.h5'
   f = h5py.File(file,'r+')
    
   x = f.get('scales/x/1.0')
   x = np.array(x)
   
   z = f.get('scales/z/1.0')
   z = np.array(z) #+ 1   # -0.5-0.5 t0 0-1

   H = 2;
   
   Nx = x.shape[0]
   Ny = z.shape[0]

   dx = x[1]-x[0]

   i = Nx//2

   f1 = 50
   f2 = 79
   
   dN = (f2-f1)*50

   bi = 0
   b = 0
   
   z_s = 1E-3*H
   z_e = 5E-1*H

   j_s = (np.abs(z-z_s)).argmin()
   j_e = (np.abs(z-z_e)).argmin()

   bx = np.zeros([Nx,Ny])
   eil2 = np.zeros([Nx])
   el2 = np.zeros([dN])
   x_min = np.zeros([dN])

   t = []

   b = np.loadtxt('BL1' + optcaseid + '.txt') #bottom wall momentum BL
   #b = np.loadtxt('BL' + optcaseid + '.txt') 
   z_H = np.loadtxt('z_H' + optcaseid + '.txt')
   print(z_H)
   print(z[j_s:j_e]/H)

   fb_i  = interpolate.interp1d(z_H,b,kind ='cubic',fill_value='extrapolate') 
   b_i = fb_i(z[j_s:j_e]/H)

   for k in range(f1, f2):
       file = path + filename + str(k) + '.h5'
       print(file)
       f = h5py.File(file,'r+')
       
       wz = f.get('tasks/wz')
       wz  = np.array(wz)

       ti = f.get('scales/sim_time')
       t  = np.concatenate((t, np.array(ti)),0)

       for p in range(0,50):
           for i,xi in enumerate(x):
               #wzi = wz[p,i,:] #bottom momentum BL
               wzi = np.flip(wz[p,i,:],0) #top momentum BL
               bx[i,:] =  2*(wzi*wzi) 
#               bi =  2*(wzi*wzi)/Nx + bi
               eil2[i] = LA.norm(b_i - bx[i,j_s:j_e])
#           bi = 0
           x_min[50*(k-f1) + p] = x[np.argmin(eil2)]
           el2[50*(k-f1) + p] = np.min(eil2)
    
   #============================================================================
   # fig1, ax1 = plt.subplots(1,figsize = [5,4])
   # ax1.plot(z_H,b)
   # ax1.plot(z[j_s:j_e]/H,b_i)
   # ax1.set_xlabel('z_H')
   # ax1.set_ylabel('b')
   # ax1.set_title(case)
   # plt.show()   
   #============================================================================
    
   fig1, ax1 = plt.subplots(1,figsize = [5,4])
   ax1.plot(t,el2)
   ax1.set_xlabel('time')
   ax1.set_ylabel('L2 norm of error')
   ax1.set_title(case)
   plt.show()  
   
   fig2, ax2 = plt.subplots(1,figsize = [5,4])
   ax2.plot(t[el2<0.7],x_min[el2<0.7],'*')
#   ax2.plot(t,x_min,'*')
   ax2.set_xlabel('time')
   ax2.set_ylabel('location of optimal match')
   ax2.set_title(case)
   plt.show() 
   
   
   print(t[el2<0.7])
   print(el2[el2<0.7])
   print(x_min[el2<0.7])    
    
   return

def find_optTBL():  #finds locations and times of optimal match
   file = path + filename + str(79) + '.h5'
   f = h5py.File(file,'r+')
    
   x = f.get('scales/x/1.0')
   x = np.array(x)
   
   z = f.get('scales/z/1.0')
   z = np.array(z) #+ 1   # -0.5-0.5 t0 0-1

   H = 2;
   
   Nx = x.shape[0]
   Ny = z.shape[0]

   dx = x[1]-x[0]

   i = Nx//2

   f1 = 50
   f2 = 79
   
   dN = (f2-f1)*50

   bi = 0
   b = 0
   
   z_s = 1E-3*H
   z_e = 5E-1*H

   j_s = (np.abs(z-z_s)).argmin()
   j_e = (np.abs(z-z_e)).argmin()

   bx = np.zeros([Nx,Ny])
   eil2 = np.zeros([Nx])
   el2 = np.zeros([dN])
   x_min = np.zeros([dN])

   t = []

   b = np.loadtxt('BL_T' + optcaseid + '.txt') #bottom wall momentum BL
   #b = np.loadtxt('BL' + optcaseid + '.txt') 
   z_H = np.loadtxt('z_H' + optcaseid + '.txt')
   print(z_H)
   print(z[j_s:j_e]/H)

   fb_i  = interpolate.interp1d(z_H,b,kind ='cubic',fill_value='extrapolate') 
   b_i = fb_i(z[j_s:j_e]/H)

   for k in range(f1, f2):
       file = path + filename + str(k) + '.h5'
       print(file)
       f = h5py.File(file,'r+')
       
       T = f.get('tasks/b')
       T  = np.array(T)

       ti = f.get('scales/sim_time')
       t  = np.concatenate((t, np.array(ti)),0)

       for p in range(0,50):
           for i,xi in enumerate(x):
               Ti = T[p,i,:] #bottom momentum BL
               #Ti = np.flip(T[p,i,:],0) #top momentum BL
               bx[i,:] =  Ti - (z-1)
#               bi =  Ti/Nx + bi
               eil2[i] = LA.norm(b_i - bx[i,j_s:j_e])
#           bi = 0
           x_min[50*(k-f1) + p] = x[np.argmin(eil2)]
           el2[50*(k-f1) + p] = np.min(eil2)
    
   #============================================================================
   # fig1, ax1 = plt.subplots(1,figsize = [5,4])
   # ax1.plot(z_H,b)
   # ax1.plot(z[j_s:j_e]/H,b_i)
   # ax1.set_xlabel('z_H')
   # ax1.set_ylabel('b')
   # ax1.set_title(case)
   # plt.show()   
   #============================================================================
    
   fig1, ax1 = plt.subplots(1,figsize = [5,4])
   ax1.plot(t,el2)
   ax1.set_xlabel('time')
   ax1.set_ylabel('L2 norm of error')
   ax1.set_title(case)
   plt.show()  
   
   fig2, ax2 = plt.subplots(1,figsize = [5,4])
   ax2.plot(t[el2<0.1],x_min[el2<0.1],'*')
#   ax2.plot(t,x_min,'*')
   ax2.set_xlabel('time')
   ax2.set_ylabel('location of optimal match')
   ax2.set_title(case)
   plt.show() 
   
   
   print(t[el2<0.1])
   print(el2[el2<0.1])
   print(x_min[el2<0.1])    
    
   return

def plot_MBL():    #plots momentum boundary layer
   file = path + filename + str(79) + '.h5'
   f = h5py.File(file,'r+')

   dt = 0.25   
   
   t1 = 100
   t2 = 125  
   
   f1 = math.floor(t1/(50*dt)) 
   f2 = math.floor(t2/(50*dt))
   
   print(f1,f2)

   dN = 1 #(f2-f1+1)*50

   x = f.get('scales/x/1.0')
   x = np.array(x)
   
   z = f.get('scales/z/1.0')
   z = np.array(z) #+ 1   # -0.5-0.5 t0 0-1

   H = 2;
   
   Nx = x.shape[0]
   Ny = z.shape[0]


   dx = x[1]-x[0]

   i = Nx//2

   #k = 12
   f1 = 50 #64 #77
   bi = 0
   b = 0
   
   z_s = 1E-3*H
   z_e = 5E-1*H

   j_s = (np.abs(z-z_s)).argmin()
   j_e = (np.abs(z-z_e)).argmin()
   print(j_s,j_e)

#   th = np.zeros((f2-f1+1)*50)
   bx = np.zeros([Nx,Ny])

   t = []
   
   jj = 23 #47 #0

   for k in range(f1, f1+1):
       file = path + filename + str(k) + '.h5'
       print(file)
       f = h5py.File(file,'r+')
       
       wz = f.get('tasks/wz')
       wz  = np.array(wz)

       ti = f.get('scales/sim_time')
       t  = np.concatenate((t, np.array(ti)),0)

       for p in range(jj,jj+1):
       
           for i,xi in enumerate(x):
               #wzi = wz[p,i,:] #bottom momentum BL
               wzi = np.flip(wz[p,i,:],0) #top momentum BL
               bx[i,:] =  2*(wzi*wzi) 
               bi =  2*(wzi*wzi)/Nx + bi

#               print(bx) 
           b = bi/dN
       
#           j_max = b[0:j_e].argmax()

#           th[(k-f1)*50 + p] = z[j_max]/H     #BL thickness
           #print(th)
           bi = 0
#   print(z[j_max]/H)

   #np.savetxt("BL"+ caseid + ".txt", b[j_s:j_e])
   #np.savetxt("z_H"+ caseid + ".txt", z[j_s:j_e]/H)

   fig1, ax1 = plt.subplots(1,figsize = [5,4])
   ax1.semilogx(z[j_s:j_e]/H,b[j_s:j_e])
   ax1.set_xlabel('z/H')
   ax1.set_ylabel('(u_x)^2 + (w_z)^2')
   ax1.set_title(case)
   plt.show()

   line = '5E6_7_opt'

   b = np.loadtxt('BL1' + line + '.txt')
   z_H = np.loadtxt('z_H' + line + '.txt')
   

   fig2, ax2 = plt.subplots(1,figsize = [5,4])
   ax2.semilogx(z_H,b, label = line, dashes = [6,2])

   b = np.loadtxt('BL' + line + '.txt')

   ax2.semilogx(z_H,b, label = '5E6_7_opt2', dashes = [2,2])
   
   ii = np.argmin(np.abs(x-15.078125))

   j = np.arange(ii,ii+1,1) 
   
   print(x[j])  

   for i,xi in enumerate(x[j]):
       ax2.semilogx(z[j_s:j_e]/H,bx[j[i],j_s:j_e], label = 'x = ' + str(x[j[i]]) + ' at t =' + str(t[jj]) )
   #ax2.set_xlabel('time')
   ax2.set_xlabel('z/H')
   ax2.set_ylabel('(u_x)^2 + (w_z)^2')
   ax2.legend(loc='upper left')
   ax2.set_title(case)
   plt.show()

   return

def plot_TBL():    #plots momentum boundary layer
   file = path + filename + str(79) + '.h5'
   f = h5py.File(file,'r+')

   dt = 0.25   
   
   t1 = 100
   t2 = 125  
   
   f1 = math.floor(t1/(50*dt)) 
   f2 = math.floor(t2/(50*dt))
   
   print(f1,f2)

   dN = 1 #(f2-f1+1)*50

   x = f.get('scales/x/1.0')
   x = np.array(x)
   
   z = f.get('scales/z/1.0')
   z = np.array(z) #+ 1   # -0.5-0.5 t0 0-1

   H = 2;
   
   Nx = x.shape[0]
   Ny = z.shape[0]


   dx = x[1]-x[0]

   i = Nx//2

   #k = 12
   f1 = 77
   bi = 0
   b = 0
   
   z_s = 1E-3*H
   z_e = 5E-1*H

   j_s = (np.abs(z-z_s)).argmin()
   j_e = (np.abs(z-z_e)).argmin()
   print(j_s,j_e)

#   th = np.zeros((f2-f1+1)*50)
   bx = np.zeros([Nx,Ny])

   t = []
   
   jj = 4 #23 #47 #0

   for k in range(f1, f1+1):
       file = path + filename + str(k) + '.h5'
       print(file)
       f = h5py.File(file,'r+')
       
       T = f.get('tasks/b')
       T  = np.array(T)

       ti = f.get('scales/sim_time')
       t  = np.concatenate((t, np.array(ti)),0)

       for p in range(jj,jj+1):
       
           for i,xi in enumerate(x):
               Ti = T[p,i,:] #bottom momentum BL
               #wzi = np.flip(wz[p,i,:],0) #top momentum BL
               bx[i,:] =  Ti - (z - 1)
               bi =  Ti/Nx + bi

#               print(bx) 
           b = bi/dN
       
#           j_max = b[0:j_e].argmax()

#           th[(k-f1)*50 + p] = z[j_max]/H     #BL thickness
           #print(th)
           bi = 0
#   print(z[j_max]/H)

   #np.savetxt("BL"+ caseid + ".txt", b[j_s:j_e])
   #np.savetxt("z_H"+ caseid + ".txt", z[j_s:j_e]/H)

   fig1, ax1 = plt.subplots(1,figsize = [5,4])
   ax1.semilogx(z[j_s:j_e]/H,b[j_s:j_e])
   ax1.set_xlabel('z/H')
   ax1.set_ylabel('T(xi,z)')
   ax1.set_title(case)
   plt.show()

   line = '5E6_7_opt'

   b = np.loadtxt('BL_T' + line + '.txt')
   z_H = np.loadtxt('z_H' + line + '.txt')
   

   fig2, ax2 = plt.subplots(1,figsize = [5,4])
   ax2.semilogx(z_H,b, label = line, dashes = [6,2])

#   b = np.loadtxt('BL' + line + '.txt')

#   ax2.semilogx(z_H,b, label = '5E6_7_opt2', dashes = [2,2])
   
   ii = np.argmin(np.abs(x-15.078125))

   j = np.arange(ii,ii+1,1) 
   
   print(x[j])  

   for i,xi in enumerate(x[j]):
       ax2.semilogx(z[j_s:j_e]/H,bx[j[i],j_s:j_e], label = 'x = ' + str(x[j[i]]) + ' at t =' + str(t[jj]) )
   #ax2.set_xlabel('time')
   ax2.set_xlabel('z/H')
   ax2.set_ylabel('(u_x)^2 + (w_z)^2')
   ax2.legend(loc='upper left')
   ax2.set_title(case)
   plt.show()

   return

def plot_energy_spectra():   # Note: computes thermal energy from perturbation temperature
   file = path + filename + str(3) + '.h5'
   f = h5py.File(file,'r+')

   dt = 0.25   
   
   t1 = 125
   t2 = 187
   
   f1 = math.floor(t1/(50*dt)) 
   f2 = math.floor(t2/(50*dt))

   print(f1,f2)

   dN = 1 #(f2-f1+1)*50

   x = f.get('scales/x/1.0')
   x = np.array(x)
   
   z = f.get('scales/z/1.0')
   z = np.array(z)
   
   Nx = x.shape[0]
   Ny = z.shape[0]

   T_k  = np.zeros([(Nx)//4,Ny], dtype=np.complex_)
   u_k  = np.zeros([(Nx)//4,Ny], dtype=np.complex_)
   w_k  = np.zeros([(Nx)//4,Ny], dtype=np.complex_)
   
   TE_k  = np.zeros([(Nx)//4], dtype=np.complex_)
   KE_k  = np.zeros([(Nx)//4], dtype=np.complex_)
   ThE_k  = np.zeros([(Nx)//4], dtype=np.complex_)
   
   TE_k = 0;
   KE_k = 0;
   ThE_k = 0;
   
   jj = Ny//2  #ymid

   z1 = -0.8
   jj = (np.abs(z-z1)).argmin()
   #jj = Ny//10  #y location at which spectra is plotted

   ThE_avg = 0;
   KE_avg = 0;
   ThE_avg_y = 0;
   KE_avg_y = 0;
   
   for k in range(f1,f1 + 1):
       file = path + filename + str(k) + '.h5'
       print(k)
       f = h5py.File(file,'r+')
   
       u = f.get('tasks/u')
       u  = np.array(u)

       w = f.get('tasks/w')
       w  = np.array(w)

       T = f.get('tasks/b')
       T  = np.array(T)

       ti = f.get('scales/sim_time')

       ThE = 0.5*np.multiply(T,T)
       KE = 0.5*(np.multiply(u,u) + np.multiply(w,w))

       ThE_avg  = ThE_avg  + ThE[:,:,:].sum(0)/dN
       KE_avg   =  KE_avg  +  KE[:,:,:].sum(0)/dN


       for ii in range(1,2):
           fft = np.fft.fft(T[ii,:,jj])/(Nx)
#           fft = fft[range(Nx//4)]
           T_k = fft
           
           fft = np.fft.fft(u[ii,:,jj])/(Nx)
#           fft = fft[range(Nx//4)]
           u_k = fft
           
           fft = np.fft.fft(w[ii,:,jj])/(Nx)
#           fft = fft[range(Nx//4)]
           w_k = fft
   
           TE_k  = TE_k + 0.5*(u_k*np.conj(u_k) + w_k*np.conj(w_k) + T_k*np.conj(T_k))/dN
           KE_k  = KE_k + 0.5*(u_k*np.conj(u_k) + w_k*np.conj(w_k))/dN
           ThE_k = ThE_k + 0.5*(T_k*np.conj(T_k))/dN
   
        
   ThE_avg_y = ThE_avg[:,jj].sum(0)/Nx
   KE_avg_y = KE_avg[:,jj].sum(0)/Nx

   k = np.arange(Nx)
   
   freq = k*alpha
#   freq = freq[range((Nx)//4)] 
   
   np.savetxt("ThE"+ caseid + ".txt", np.real(ThE_k))
   np.savetxt("TE"+ caseid + ".txt", np.real(TE_k))
   np.savetxt("KE"+ caseid + ".txt", np.real(KE_k))

   np.savetxt("kx"+ caseid + ".txt", freq)

   print(ThE_avg_y)
   print(ThE.shape[0])
   
   fig2, ax2 = plt.subplots(1,figsize=[7,4])
   ax2.loglog(freq,np.real(TE_k),'r', label = 'Total Energy')
   ax2.loglog(freq,np.real(KE_k), 'b',label = 'Kinetic Energy')
   ax2.loglog(freq,np.real(ThE_k), 'g', label = 'Thermal Energy')
   ax2.set_xlabel('k')
   ax2.set_ylabel('E_k')
   ax2.set_title(case)
#   ax2.set_ylim((1E-9,1E-1))
   #ax2.legend(bbox_to_anchor=(0.1, 0.), loc='upper right')
   ax2.legend(loc='upper right')
  
   plt.tight_layout()
   plt.show()

   X,Z = np.meshgrid(x,z)
   N = 50

#   fig, axarr = plt.subplots(2,2, figsize=(12,7))
#   levels = np.linspace(0,0.7,21)
#   levels1 = np.linspace(0,0.2,21)
#   axarr[0,0].contourf(X, Z, ThE[N-20,:,:].T,levels,cmap=cm.seismic)
#   axarr[0,0].set_title('Thermal Energy at t = %s'%(ti[N-20]))
#   axarr[1,0].contourf(X, Z, ThE_avg.T,levels,cmap=cm.seismic)
#   axarr[1,0].set_title('Averaged Thermal Energy')
#   axarr[0,1].contourf(X, Z, KE[N-20,:,:].T,levels1)
#   axarr[0,1].set_title('Kinetic Energy at t = %s'%(ti[N-20]))
#   axarr[1,1].contourf(X, Z, KE_avg.T,levels1)
#   axarr[1,1].set_title('Averaged Kinetic Energy')

   plt.tight_layout()
   plt.show()
   
   return 

caseid = '5E6_7_20'
optcaseid = '5E6_7_opt'
path = '/Volumes/Chest/snapshots_' + caseid + '/'
#path = '../snapshots_1E7_10_32/'
#path = '../1E7_4_opt12.01_wom/'
filename = 'snapshots_'+ caseid+ '_s'
case = 'Ra = 5E6, Pr = 7, Box size = 20'
#alpha = 0.4714574564629509E+001 # Wavenumber we're working with
#alpha = 0.1560021496301813E+002 # Wavenumber we're working with
alpha = 0.5080785292900329E+001; #5E6, 7
 
#onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
#onlyfiles = np.array(onlyfiles)
#nfiles = onlyfiles.shape[0]

#print(nfiles)
#spatial_correlation()
#opt_comparison()
#find_optTBL()
#plot_TBL()
#plot_energy_spectra()
#plot_Nu()
track_localNu(-2.421875)
size = [32]
#size = [64, 32, 24, 12, 6, 5, 3]
E_type = 'ThE'
#plot_energyspectra()
E_type = 'KE'
#plot_energyspectra()
E_type = 'TE'
#plot_energyspectra()

#plot_contours()
#Nu_avg = [14.583494728539879 ,14.38040974149913, 14.412647083759014, 14.781802371717072, 15.302233353679274, 15.634397375288362 ]
#Nu_avg = [14.83824719915766, 14.576552724393101, 14.93889799347675, 16.220887956354023, 16.71374860004233, 16.22138431031129  ]
#fig1, ax1 = plt.subplots(1,figsize = [5,4])
#ax1.plot(size,Nu_avg)
#ax1.set_xlabel('size')
#ax1.set_ylabel('Nu_avg')
#plt.show()
