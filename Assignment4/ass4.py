"""
Assignment 4, EE2703 - Applied Lab Programming

Dhruv Goyal, EE19B077

Usage example: python3 EE2703_ASSIGN4_EE19B077.py
"""

# importing libraries
import scipy.integrate as integrate
import math
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

#defining functions

#exponential
def exp(x):
	y = np.exp(x)
	return y

#Cos(Cos(x))
def coscos(x):
	y = np.cos(np.cos(x))
	return y

#Function to find out the fourier coeffcients for exp(x)
def expcoeff():
	i = 1
	j = 1
	n = 25

	def u(x,k):
		return np.exp(x)*np.cos(k*x)		#defining u as exp(x)*Cos(kx) to pass through the quad function with arguments as k

	def v(x,k):
		return np.exp(x)*np.sin(k*x)		#defining v as exp(x)*Sin(kx) to pass through the quad function with arguments as k


	exp_coeff[0] = integrate.quad(u,0,2*pi,args=(0))[0]		#Calculating A0 using u(x,k)
	exp_coeff[0] *= (1/(2*pi))

	while (i<=2*n):
		exp_coeff[i] = integrate.quad(u,0,2*pi,args=(j))[0]		#Calculating the rest of the Ak coeffcients
		exp_coeff[i] *= (1/pi)
		i = i+1
		exp_coeff[i] = integrate.quad(v,0,2*pi,args=(j))[0]		#Calculating the rest of the Bk coeffcients
		exp_coeff[i] *= (1/pi)
		j = j+1
		i = i+1


#Function to find out the fourier coeffcients for Cos(Cos(x)) with procedure the same as for exp(x)
def coscoeff():
	i = 1
	j = 1
	n = 25

	def u1(x,k):
		return coscos(x)*np.cos(k*x)

	def v1(x,k):
		return coscos(x)*np.sin(k*x)

	cos_coeff[0] = integrate.quad(u1,0,2*pi,args=(0))[0]
	cos_coeff[0] *= (1/(2*pi))

	while (i<=2*n):
		cos_coeff[i] = integrate.quad(u1,0,2*pi,args=(j))[0]
		cos_coeff[i] *= (1/pi)
		i = i+1
		cos_coeff[i] = integrate.quad(v1,0,2*pi,args=(j))[0]
		cos_coeff[i] *= (1/pi)
		j = j+1
		i = i+1


#Function to plot the fourier coeffcients of exp(x) and takes input as the best fit vector found out using least sqaure method
def fourierplotexp(c0):

	plt.figure(3)																#plotting semilog graph
	plt.semilogy(arr, abs(exp_coeff),'or',label = 'integrate_coeff')			#plotting coefficents obtained using direct integration methods as red dots
	plt.semilogy(arr, abs(c0),'og',label = 'estimated_coeff', markersize=5)		#plotting coefficents obtained using estimation by least square method as green dots

	plt.title("Semilogy plot for coeff of exp(x)")								#details and labels for the graph
	plt.xlabel('n')
	plt.ylabel('An or Bn')
	plt.grid(True)
	plt.legend()
	plt.show()


	plt.figure(4)																#plotting loglog graph
	plt.loglog(arr, abs(exp_coeff),'or',label = 'integrate_coeff')				#plotting coefficents obtained using direct integration methods as red dots
	plt.loglog(arr, abs(c0),'og',label = 'estimated_coeff', markersize=5)		#plotting coefficents obtained using estimation by least square method as green dots

	plt.title("loglog plot for coeff of exp(x)")								#details and labels for the graph
	plt.xlabel('n')
	plt.ylabel('An or Bn')
	plt.grid(True)
	plt.legend()
	plt.show()


#Function to plot the fourier coeffcients of Cos(Cos(x)) and takes input as the best fit vector found out using least sqaure method
def fourierplotcos(c1):
	plt.figure(5)																#plotting semilog graph
	plt.semilogy(arr,abs(cos_coeff),'or',label = 'integrate_coeff')				#plotting coefficents obtained using direct integration methods as red dots
	plt.semilogy(arr, abs(c1),'og',label = 'estimated_coeff',markersize=5)		#plotting coefficents obtained using estimation by least square method as green dots


	plt.title("Semilogy plot for coeff of coscos(x)")							#details and labels for the graph
	plt.xlabel('n')
	plt.ylabel('An or Bn')
	plt.grid(True)
	plt.legend()	
	plt.show()
	
	plt.figure(6)																#plotting loglog graph
	plt.loglog(arr,abs(cos_coeff),'or',label = 'integrate_coeff')				#plotting coefficents obtained using direct integration methods as red dots
	plt.loglog(arr, abs(c1),'og',label = 'estimated_coeff',markersize=5)		#plotting coefficents obtained using estimation by least square method as green dots

	plt.title("loglog plot for coeff of coscos(x)")								#details and labels for the graph
	plt.xlabel('n')
	plt.ylabel('An or Bn')
	plt.grid(True)
	plt.legend()	
	plt.show()

#To find the fourier series coefficients of exp(x) using least aquare approach, finding the error and plot the coefficients and the resulting function along with the original one
def leastsquareexp(x,A):
	b = exp(x)
	c0 = np.linalg.lstsq(A, b, rcond=None)[0]      # the ’[0]’ is to pull out the# best fit vector. lstsq returns a list.

	err = np.subtract(exp_coeff,c0)      #To find deviation between integrated and estimated coeffcients
	max_exp=max(abs(err))				 #Finding maximum absolute error

	print(err)								#printing out the error matrix
	print("max error exp(x) =", end = " ")
	print(max_exp)							#printing out the maximum error found
	

	plt.figure(1)							#Plotting semilog plot of the original function from _2*pi to 4*pi and the estimated function using fourier coeffecients from 0 to 2*pi
	plt.grid(True)
	plt.semilogy(X,exp(X), "r", label = 'Original Func')
	plt.semilogy(x,np.dot(A,c0),'og',label = 'Estimated Func', markersize=2)
	plt.xlabel(r'x$\rightarrow$',fontsize=12)
	plt.ylabel(r'exp(x)$\rightarrow$',fontsize=12)
	plt.legend(loc='upper right')
	plt.show()

	return c0								#returning the best fit vector

#To find the fourier series coefficients of Cos(Cos(x)) using least aquare approach, finding the error and plot the coefficients and the resulting function along with the original one
def leastsquarecos(x,A):
	b = coscos(x)
	c1=np.linalg.lstsq(A,b, rcond=None)[0]      # the ’[0]’ is to pull out the# best fit vector. lstsq returns a list.

	err = np.subtract(cos_coeff,c1)      		#To find deviation between integrated and estimated coeffcients
	max_exp=max(abs(err))						#Finding maximum absolute error
	
	print(err)									#printing out the error matrix
	print("max error cos(cos(x) =", end = " ")
	print(max_exp)								#Maximum deviation


	plt.figure(2)								#Plotting semilog plot of the original function from _2*pi to 4*pi and the estimated function using fourier coeffecients from 0 to 2*pi
	plt.grid(True)
	plt.semilogy(X,coscos(X),"r", label = 'Original Func')
	plt.semilogy(x,np.dot(A,c1),'og',label = 'Estimated Func', markersize=2)
	plt.xlabel(r'x$\rightarrow$',fontsize=12)
	plt.ylabel(r'cos(cos(x)$\rightarrow$',fontsize=12)
	plt.legend(loc='upper right')
	plt.show()

	return c1									#returning the best fit vector


"""DRIVER COMMANDS"""

pi = math.pi                                 #defining pi
exp_coeff = np.zeros(51,dtype = float)       #Storing oefficients of exp(x) using integration method
cos_coeff = np.zeros(51,dtype = float)       #Storing coefficients of cos(cos(x)) using integration method
arr = np.arange(51)

X = np.linspace((-2*pi),4*pi,400)			#defines X from -2*pi to 4*pi in 400 steps
x=np.linspace(0,2*pi,401)					#defines X from 0 to 2*pi in 401 steps
x=x[:-1] # drop last term to have a proper periodic integral

A=np.zeros((400,51))     # allocate space for A
A[:,0]=1              # col 1 is all ones
for k in range(1,26):
	A[:,2*k-1]=np.cos(k*x) # cos(kx) column
	A[:,2*k]=np.sin(k*x)   # sin(kx) column


#Calling all the defined functions
expcoeff()
coscoeff()
c0 = leastsquareexp(x,A)
c1 = leastsquarecos(x,A)
fourierplotexp(c0)
fourierplotcos(c1)
	



