from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import numpy.linalg as linalg 
import math


def calc(l):
	R = r_vector - r_prime[l,:]
	R_mod = linalg.norm(R, axis=1)

	R_mod = transpose(R_mod, axes=[1,0,2])

	k = 0.1
	F = (cos(phi[l])/R_mod)*(e**(k*R_mod*-1j)).reshape(1,9000)

	A = np.transpose(F)*dl[l,:].reshape(1,3)
	return A


def estimateExponent(vec):
	G = np.array([[1, log(x)] for x in range(1,len(vec)+1)])	#Defining array with ones and x for both the cases
	x = linalg.lstsq(G, np.log(vec), rcond = None)[0]			#using least square approach to find values of logA and B, [0] gives best fit
	return x[0], x[1]											#returning these values
	

x = np.linspace(-1,1,3)
y = np.linspace(-1,1,3)
z = np.linspace(1,1000,1000)

u0 = 4*pi*(1e-7)

xx,yy,zz = np.meshgrid(x,y,z)

radius = 10
L = 100


phi = np.linspace(0,2*pi,L+1)
phi = phi[:-1]

r_prime = np.c_[radius*cos(phi), radius*sin(phi), np.zeros(100)]

r_vector = np.c_[xx.reshape(-1), yy.reshape(-1), zz.reshape(-1)]

dl = np.c_[2*pi*(-sin(phi)*cos(phi))/10, 2*pi*(sin(phi)*cos(phi))/10, np.zeros(100)]

I_vector = np.c_[4*pi*cos(phi)*(-sin(phi))/u0, 4*pi*(cos(phi))**2/u0, np.zeros(100)]

quiver(r_prime[:,0], r_prime[:,1], I_vector[:,0], I_vector[:,1])
legend(loc = "best")
grid(True)
show()

A = 0

for l in range(100):
	A+= calc(l)

A_x = A[:,0]
A_y = A[:,1]

A_x = A_x.reshape(3,3,1000)
A_y = A_y.reshape(3,3,1000)

B=(A_y[2,1,:]-A_x[1,2,:]-A_y[0,1,:]+A_x[1,0,:])/(4)

c,b = estimateExponent(B)

estimate = (exp(c))*(z**b)

loglog(z,abs(B))
loglog(z,estimate)
grid(True)
show()



