"""
Assignment 5, EE2703 - Applied Lab Programming

Dhruv Goyal, EE19B077

Usage example: python3 EE2703_ASSIGN5_EE19B077.py
"""

# importing libraries
import mpl_toolkits.mplot3d.axes3d as p3 #module used for surface plot
from sys import argv
import numpy as np
import scipy.linalg as lin
from pylab import *
import matplotlib.pyplot as plt

# Defining Functions
def plot_initial_contour():
	plt.figure(1)
	plt.contourf(X,Y, phi)														#Contourfill plot
	plt.scatter(ii[0] - (Nx - 1) / 2, ii[1] - (Ny - 1) / 2, color = 'r')		#Red dots to show the region where the lead wire is
	plt.colorbar()																#Showing colourbar
	plt.title('Initialised Contour Plot')										#Giving x and y labels as well as the title of plot
	plt.xlabel(r'x$\rightarrow$',fontsize=12)
	plt.ylabel(r'y$\rightarrow$',fontsize=12)
	plt.show()


# Applying Boundary Condition
def bound(phi,ii):
    phi[:,0]=phi[:,1]                    # Left Boundary
    phi[:,Nx-1]=phi[:,Nx-2]              # Right Boundary
    phi[0,:]=phi[1,:]                    # Top Boundary
    phi[Ny-1,:]=0
    phi[ii]=1.0
    return phi

#Finding error and updating phi
def err_each_iteration(Niter,phi,err,ii):
    for k in range(Niter):
        old_phi = phi.copy()  #Making copy of the phi array
        phi[1:-1,1:-1] = 0.25*(old_phi[1:-1,0:-2]+ old_phi[1:-1,2:]+ old_phi[0:-2,1:-1] + old_phi[2:,1:-1]) 	#Updating phi as average of adjacent values
        phi = bound(phi,ii)																					 	#Applying boundary conditions
        err[k] = np.max(np.abs(phi-old_phi))																	#Updating error array

#Plotting the Error
def err_plotting(Niter,err):
    #plotting Error and deining labels
    plt.figure(2)
    plt.title('Error vs iterations')
    plt.xlabel(r'No of Iterations $\rightarrow$')
    plt.ylabel(r'Error$\rightarrow$')
    plt.plot(range(Niter),err,'-r',markersize=3,label='original')
    plt.legend()
    plt.show()

    #plotting Error on semilog
    plt.figure(3)
    plt.title('Error vs iterations on semi log scale')
    plt.xlabel(r'No of Iterations $\rightarrow$')
    plt.ylabel(r'Error$\rightarrow$')
    plt.semilogy(range(Niter),err,'-r',markersize=3,label='original')
    plt.legend()
    plt.show()

    #plotting Error on loglog
    plt.figure(4)
    plt.title('Error vs iterations on loglog scale')
    plt.xlabel(r'No of Iterations $\rightarrow$')
    plt.ylabel(r'Error$\rightarrow$')
    plt.loglog((np.asarray(range(Niter))+1),err)
    plt.loglog((np.asarray(range(Niter))+1)[::50],err[::50],'ro')
    plt.legend(["real","every 50th value"])
    plt.show()

#Estimation for the exponential
def estimateExponent(n,vec):
	A = np.array([[1, x] for x in range(500*n,len(vec)+500*n)])	#Defining array with ones and x for both the cases
	x = lin.lstsq(A, np.log(vec))[0]							#using least square approach to find values of logA and B, [0] gives best fit
	return x[0], x[1]											#returning these values
	
#Plotting best fit and error
def plot_error_fit():
	#semilog scale
	plt.figure(5)
	plt.title("Best fit for error on semilog scale")
	plt.xlabel(r"Iterations$\longrightarrow$")
	plt.ylabel(r"Error$\longrightarrow$")
	x = np.asarray(range(Niter))+1
	plt.semilogy(x,err,label="Error")
	plt.semilogy(x[::50],np.exp(A1+B1*np.asarray(range(Niter)))[::50],'ro',label="fit1")      #plotting in jumps of 50 as is asked
	plt.semilogy(x[::50],np.exp(A2+B2*np.asarray(range(Niter)))[::50],'go',markersize = 4,label="fit2")
	plt.legend()
	plt.show()

	#loglog scale
	plt.figure(6)
	plt.title("Best fit for error on loglog scale")
	plt.xlabel(r"Iterations$\longrightarrow$")
	plt.ylabel(r"Error$\longrightarrow$")
	x = np.asarray(range(Niter))+1
	plt.loglog(x,err,label="Error")
	plt.loglog(x[::50],np.exp(A1+B1*np.asarray(range(Niter)))[::50],'ro',label="fit1")      #plotting in jumps of 50 as is asked
	plt.loglog(x[::50],np.exp(A2+B2*np.asarray(range(Niter)))[::50],'go',markersize = 4,label="fit2")
	plt.legend()
	plt.show()

#Plotting the Contour Plot of Potential
def Contour_plotting(X,Y,Nx,Ny):
    #plotting 2d contour of final potential
    fig1=plt.figure(7)  
    plt.title(r'2D Contour plot of potential')
    plt.plot(ii[1]-(Nx-1)/2,ii[0]-(Ny-1)/2,'ro')					#plotting red dots where the lead wire is present/unit potential 
    plt.contourf(Y,X[::-1],phi)
    plt.ylabel(r'y$\rightarrow$')
    plt.xlabel(r'x$\rightarrow$')
    plt.colorbar()
    plt.show()

#Plotting 3D Plot of Potential
def poten3D_plotting(X,Y,phi):
    #plotting 3d contour of final potential
    fig1=plt.figure(8)     															# open a new figure
    ax=p3.Axes3D(fig1) 																# Axes3D is the means to do a surface plot
    surf = ax.plot_surface(Y, X, phi.T, rstride=1, cstride=1, cmap=plt.cm.jet)
    ax.set_title(r'The 3-D surface plot of the potential')
    plt.xlabel(r'x$\rightarrow$',fontsize=15)
    plt.ylabel(r'y$\rightarrow$',fontsize=15)
    ax.set_zlabel(r'$\phi\rightarrow$',fontsize=15)
    plt.show()

#PLotting the current direction
def current_plotting(X,Y,Jx,Jy,Nx,Ny):
    #plotting current density
    plt.figure(9)
    plt.title("Vector plot of current flow")
    plt.quiver(-x,y,-Jx[::-1,:],-Jy[::-1,:],scale=4) 	#quiver helps plotting arrows in direction of vector
    plt.plot(ii[1]-(Nx-1)/2,ii[0]-(Ny-1)/2,'ro') 		#plotting red dots for unit potential
    plt.title('contour of current density')
    plt.xlabel(r'x$\rightarrow$',fontsize=15)
    plt.ylabel(r'y$\rightarrow$',fontsize=15)
    plt.show()



#Driver Commands

#if arguments are given, checking if correct number of argumnets given
try:
	len(argv)==1 or len(argv)==5
except IOError:
	exit(0)


#get user inputs
if(len(argv)==5):
    Nx=int(argv[1])
    Ny=int(argv[2])
    radius=int(argv[3])
    Niter=int(argv[4])
    print("Using user provided parameters")

else:
    Nx=25                                            # size along x
    Ny=25                                            # size along y
    radius=8                                         #radius of central lead
    Niter=1500                                       #number of iterations to perform
    print("Using defaults. If u want to use your own, put all 4 parameters as commandline arguments")


#Initialising position arrays x and y
x = np.linspace(-(Nx - 1) / 2, (Nx - 1) / 2, Nx, dtype = float) 
y = np.linspace(-(Ny - 1) / 2, (Ny - 1) / 2, Ny, dtype = float)

phi = np.zeros((Nx,Ny),dtype=float) 			#creating a matrix array called phi
Y, X = np.meshgrid(y, x)						#Creating meshgrid using x and y arrays
ii = np.where(X * X + Y * Y <= radius ** 2)		#Finding indices of points within the radius
phi[ii] = 1.0									#updating these points to unit potential

#Plotting initalised potential
plot_initial_contour()

#Plotting error
err = np.zeros(Niter, dtype = float)			#defining error matrix
err_each_iteration(Niter,phi,err,ii)			#using algorithm to find potential and the error
err_plotting(Niter,err)							#plotting the error

A1, B1 = estimateExponent(0,err)				#Finding best fit using lstsq for entire error vector
A2, B2 = estimateExponent(1,err[500:])			#Finding best fit using lstsq from the 500th iteration in error vector

plot_error_fit()								#Plotting the best fit along with original curve

Contour_plotting(X,Y,Nx,Ny)						#Plotting 2-D contour plot of final potential
poten3D_plotting(X,Y,phi)						#Plotting 3-D contor plot of final potential

#initialising jx and jy arrays
Jx = np.zeros((Nx,Ny))
Jy = np.zeros((Nx,Ny))

#Finding current density using given relation
Jx[1:-1,1:-1] = 0.5*(phi[1:-1,0:-2]-phi[1:-1,2:])
Jy[1:-1,1:-1] = 0.5*(phi[0:-2,1:-1]-phi[2:,1:-1])

current_plotting(x,y,Jx,Jy,Nx,Ny)				#Plotting the vector current density
