"""
Endsemester Examination, EE2703 - Applied Lab Programming

Dhruv Goyal, EE19B077

Usage example: python3 EE2703_Endsem_E19B077.py               

In this assignment, we analyze the magnetic field along z-axis for a given distribution of current and find its decay rate using least square approach.     
"""
# Importing libraries 
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import numpy.linalg as linalg 
from scipy.linalg import lstsq
import math

#Defining constants
u0 = 4*pi*(1e-7)


# Defining functions 

# calc function calculates the value of |R| and uses it to find the value of magnetic vector potential

def calc(l):													# Takes input as an integer l to mark a given point on loop  
	R = r_vector - r_prime[l]									# R is the vector from the point on loop to the point at which A is found
	R_mod = linalg.norm(R, axis=3)								# Norm outputs the magnitude of R

	k = 0.1
	alpha = ((cos(phi[l]))/R_mod)*(e**(k*R_mod*(-1j)))

	A_x = alpha*dl[l,0]											# Finding x component of A vector
	A_y = alpha*dl[l,1]											# Finding y component of A vector

	return A_x,A_y												# Returning values


# Estimate function takes a vector as input and fits an exponential to it using lstsq
def estimateExponent(vec):
	T = np.array([[1, log(x)] for x in range(1,len(vec)+1)])	# Defining array with ones and log(x)
	x = linalg.lstsq(T, np.log(vec), rcond = None)[0]			# using least square approach to find values of logA and B, [0] gives best fit
	return x[0], x[1]											# returning values


# Function to divide the vector space into a meshgrid								
def mesh():
	x = np.linspace(-1,1,3)										# Defining x,y and z arrays using linspace
	y = np.linspace(-1,1,3)
	z = np.linspace(1,1000,1000)
	xx, yy, zz = np.meshgrid(x, y, z, indexing = 'ij')			# Indexing as 'ij' instead of default which is 'ji'
	return x, y, z, xx, yy, zz 


#Plotting scatter plot of points on loop segments
def plot_scatter():												
	scatter(r_prime[:,0],r_prime[:,1], color = "red", label="Segment Points")
	title("Scatter plot of points on segments of the loop")
	ylabel(r'y$\rightarrow$')
	xlabel(r'x$\rightarrow$')
	legend(loc="best")
	grid(True)
	show()

#Plotting quiver plot of current distribution
def plot_quiver():												
	quiver(r_prime[:,0], r_prime[:,1], I_vector[:,0], I_vector[:,1],color = "blue",label="Current element")
	title("Vector plot of the current distribution of the loop")
	ylabel(r'y$\rightarrow$')
	xlabel(r'x$\rightarrow$')
	legend(loc = "best")
	grid(True)
	show()



# Driver Commands

x,y,z,xx,yy,zz = mesh()											# Using meshgrid to define a vector space

r_vector = np.zeros((3,3,1000,3))								# Using a 4-D array to define a (3,3,1000) 3-D array of x,y,z values
r_vector[:,:,:,0] = xx
r_vector[:,:,:,1] = yy
r_vector[:,:,:,2] = zz


radius = 10 													# Value of radius of loop
N = 100 														# No. of sections of loop


phi = np.linspace(0,2*pi,N+1)									# Angle phi broken into N+1 parts
phi = (phi[:-1]+phi[1:])/2										# The points are set to the mid points of line segments


r_prime = np.c_[radius*cos(phi), radius*sin(phi), np.zeros(100)]							# Position vector array of points in the loop	

dl = np.c_[2*pi*(-sin(phi))/10, 2*pi*(cos(phi))/10, np.zeros(100)]							# Defining element vector along the loop

I_vector = np.c_[4*pi*cos(phi)*(-sin(phi))/u0, 4*pi*(cos(phi))**2/u0, np.zeros(100)]		# Defining the value of current element for dl


plot_scatter()													#Plotting the loop segments using scatter plot
plot_quiver()													#Plotting the current element vectors using quiver plot


A_x, A_y = 0,0 													# Initializing A_x, A_y to 0

for l in range(0,N):											# Using for loop to evaluate and sum the values of A_x and A_y using calc
	A_x += calc(l)[0]
	A_y += calc(l)[1]

B=(A_y[1,0,:]-A_x[0,1,:]-A_y[-1,0,:]+A_x[0,-1,:])/4				# Evaluating the magnetic field (off-axis)

c,b = estimateExponent(abs(B))									# Using least square approach to find c,b in c*z^b

print("exp(c) = {}, b = {}".format(exp(c), b))					# Printing out the values

best_fit = (exp(c))*(z**b)										# Defining the best fit function 


#Plotting magnetic field vector on a log-log plot
title(r'Log Log plot of the $|\vec{B_z}|$ vs z(Off-axis)')
loglog(z,abs(B),'b', label = 'Magnetic field')                
legend(loc = "best")
grid(True)
ylabel(r'$log(|\vec{B_z}|)\rightarrow$')
xlabel(r'log(z)$\rightarrow$')
show()

#Plotting magnetic field vector along with best fit as (c*z^b) on a log-log plot 
plt.title(r'Log Log plot of the $|\vec{B_z}|$ vs z(Off-axis) with best fit')
plt.loglog(z,abs(B),'b', label = 'Original')                                      
plt.loglog(z,best_fit,'r',label = 'Least square fit')                                     
plt.legend(loc = "best")
plt.grid(True)
plt.ylabel(r'$log(|\vec{B_z}|)\rightarrow$')
plt.xlabel(r'log(z)$\rightarrow$')
plt.show()