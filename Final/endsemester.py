"""
Program to Find the Magnetic field along z-axis for a given current value

Submitted by: EE19B133 Puneet Sangal
Written on: 23/05/2021, Sunday

Usage example: python EE19B133_ENDSEM_EE2703.py
Output: Gives different graph for current and magnetic field varying along z-axis.
"""

# importing the necessary files for plotting and math computing
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as linalg
import math

#defining the calc function which gives absolute value of R for a given point on loop and A vector
def calc(l):
	R = r_vec - r_pr[l]                                                               #finding R vector
	R_mod = linalg.norm(R, axis=-1)                                                   #finding modulus of R vector
	R_mod = np.transpose(R_mod, axes=[1,0,2])

	k = 0.1
	F = (abs(np.cos(phi[l]))/R_mod)*(np.exp(k*R_mod*(-1j)))                           #Finding the terms of A

	A_x = F*dl[l,0]
	A_y = F*dl[l,1]

	return A_x,A_y                                                                    #return A_x and A_y

#Finding C and B for best fit
def estimate(vec):
	G = np.array([[1, np.log(x)] for x in range(1,len(vec)+1)])
	x = linalg.lstsq(G, np.log(vec), rcond = None)[0]                                 #Using numpy.linalg.lstsq for estimate
	return x[0], x[1]


# Main Function

# Defining the points taken in 3-D plot
x = np.linspace(-1,1,3)
y = np.linspace(-1,1,3)
z = np.linspace(1,1000,1000)

# Defining the Value of magnetic permeability
u0 = 4*np.pi*(1e-7)

#Getting all points in 3-D Plane as 3-D array as [3,3,1000]
X,Y,Z = np.meshgrid(x,y,z)

#Defining Vector for each point in 3-D array
r_vec = np.zeros((3,3,1000,3))
r_vec[:,:,:,0] = X
r_vec[:,:,:,1] = Y
r_vec[:,:,:,2] = Z

#Defining Parameters of loop
radius = 10
L = 100

#Defining value of angle of 100 points taken on loop
phi = np.linspace(0,2*np.pi,L+1)
phi = phi[:-1]

#Defining vector for each point taken on the loop
r_pr = np.c_[radius*np.cos(phi), radius*np.sin(phi), np.zeros(100)]

#Defining the current element vector for each point on loop
dl = np.c_[2*np.pi*(-np.sin(phi))/10, 2*np.pi*(np.cos(phi))/10, np.zeros(100)]

#Defining the Current vector along current element vector
I_vec = np.c_[4*np.pi*np.cos(phi)*(-np.sin(phi))/u0, 4*np.pi*(np.cos(phi))**2/u0, np.zeros(100)]

# Plotting the Loop taken in X-Y plane along with current vector at each point
dx, dy = np.gradient(I_vec)
color = np.sqrt((abs(dx[:,0]+2)/2)*2 + (abs(dy[:,1]+2)/2)*2)                                                           #defining color for each vector
plt.title("Vector plot of current flow on the loop")
plt.quiver(r_pr[:,0], r_pr[:,1], I_vec[:,0], I_vec[:,1],color , label = "current vector",color = 'yellow')             #doing quiver plot
plt.legend(loc = "best")
plt.grid(True)
plt.ylabel(r'y$\rightarrow$')
plt.xlabel(r'x$\rightarrow$')
plt.show()

A_x, A_y = 0,0

#Finding the values of A_x and A_y using Calc() function
for l in range(100):
	A_x += calc(l)[0]
	A_y += calc(l)[1]

#Finding Magnetic field along z-Axis
B=(A_y[2,1,:]-A_x[1,2,:]-A_y[0,1,:]+A_x[1,0,:])/4

# Plotting The Magnetic Field on log-log scale
plt.title(r'Log Log plot of the $|\vec{B_z}|$ vs z')
plt.loglog(z,abs(B),'b', label = 'original')                                       # Plotting log-log plot
plt.legend(loc = "best")
plt.grid(True)
plt.ylabel(r'$|\vec{B_z}|\rightarrow$')
plt.xlabel(r'z$\rightarrow$')
plt.show()

#Finding estimate of C and B for B = c*z^b
c,b = estimate(B)
print(np.exp(c),b)
estimate = (np.exp(c))*(z**b)

# Plotting The Magnetic Field on log-log scale along with estimate graph
plt.title(r'Log Log plot of the $|\vec{B_z}|$ vs z')
plt.loglog(z,abs(B),'b', label = 'original')                                      #Plotting Log-log plot
plt.loglog(z,estimate,'r',label = 'estimate')                                     #Plotting Log-Log plot on an estimate
plt.legend(loc = "best")
plt.grid(True)
plt.ylabel(r'$|\vec{B_z}|\rightarrow$')
plt.xlabel(r'z$\rightarrow$')
plt.show()
