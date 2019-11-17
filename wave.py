#!/usr/bin/env python


import matplotlib.pyplot as plt
import numpy as np 
import os
import subprocess
import glob
#FINITE DIFFERENCE 1D WAVE EQUATION

#########################################
#define model parameters

l=10000 #length [m]
nsteps=2000 #number of timesteps
dx=25 #[m]
dt=0.005 # [sec]
elements=int(l/dx)
density=2800 #[kg/m^3]
velocity=3000 #[m/s]

####################################################
#setting the arrays

U=np.ndarray(elements,dtype=float);U[:]=0
Upast=np.ndarray(elements,dtype=float);Upast[:]=0
Ufuture=np.ndarray(elements,dtype=float);Ufuture[:]=0
lam=np.ndarray(elements,dtype=float);
rho=np.ndarray(elements,dtype=float); rho[:]=density
vel=np.ndarray(elements,dtype=float); vel[:]=velocity
lam[:]=rho[:]*vel[:]**2
epsilon=np.ndarray(elements,dtype=float)
sig=np.ndarray(elements,dtype=float)
RHS=np.ndarray(elements,dtype=float)
x=np.ndarray(elements,dtype=float)

################################################################
# NOTE: this model is homogeneous. Interfaces ,i.e., changing medium
# needs to be set menualy

#################################################################
#initial conditions
a=5.5e-6
x0=2000
for i in range(0,elements):
	x[i]=(i-1)*dx
	U[i]=np.exp(-a*(x[i]-x0)**2)
	Upast[i]=np.exp(-a*(x[i]-(x0-vel[i]*dt))**2)


###################################################################
#INTEGRATION SPACE AND TIME
for i in range(0,nsteps):
	for j in range(0,elements-1):
		epsilon[j]=(U[j+1]-U[j])/dx
		sig[j]=lam[j]*epsilon[j]

	for j in range(1,elements-1):	
		RHS[j]=(sig[j]-sig[j-1])/dx

	Ufuture[:]=(RHS[:]/rho[:])*dt**2+2*U[:]-Upast[:]
	Ufuture[0]=0
	Ufuture[elements-1]=0

	#################################################
	#printing figures
	#change temp value in order to change the number of printed figures
	temp=10
	if i%temp==0:
		p_num=i/temp
		plt.figure()
		plt.plot(x,Ufuture[:]); plt.title("1D wave"); plt.xlabel('Distance [m]') 
		exec "plt.savefig('1D-wave%d.png')" %(p_num)
	#################################################



	Upast[:]=U[:]
	U[:]=Ufuture[:]

#######################################################################
#In order to make a movie use:

#movie subplots
subprocess.call(['ffmpeg' ,'-f', 'image2' ,'-r' ,'2' ,'-i' ,'1D-wave%d.png' ,'-vcodec' ,'mpeg4' ,'-y','movie.mp4'])

for f in glob.glob("1D-wave*.png"):
    os.remove(f)


