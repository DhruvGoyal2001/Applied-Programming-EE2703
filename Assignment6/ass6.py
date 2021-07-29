"""
Assignment 6, EE2703 - Applied Lab Programming

Dhruv Goyal, EE19B077

Usage example: python3 EE2703_ASSIGN6_EE19B077.py n=value M=value nk=Value u0=Value p=Value Msig=Value 
where all the arguments are optional and value is some number                

Inputs: M (no. of electrons injected) , n (len of tubelight) , nk (time steps to be taken) , u0 (velocity after which collision is possible)
        p (probability that an electron at velocity higher than u0 collides) , Msig (value of standard deviation for normal distribution of number of electrons injected)
        
Outputs: No. of electrons at a position vs position graph
         Intensity of light at a position vs position graph
         Velocity of electron at a position vs position graph
         Tabular data of Intensity vs position         
"""

#Importing libraries
from sys import argv
import numpy as np    
import matplotlib.pyplot as plt
from pylab import *
import pandas as pd    # pandas for showing in tabular form


#Defining functions

def electron_info():

	for i in range(1,nk):								

		if i==1:											# Finding indices where electron is present
			ii=np.where(xx>0)[0]							# first element of tuple returned by np.where gives the indices based on given condition

		dx[ii] = u[ii] + 0.5								# Updating values, acceleration = 1 unit
		xx[ii] = xx[ii] + dx[ii]
		u[ii] = u[ii] + 1.0

		npos = np.where(xx>=n)[0]							# Finding indices of electrons which have reached anode
		xx[npos] = 0.0										# Resetting all values for those electrons
		dx[npos] = 0.0
		u[npos] = 0.0

		kk=np.where(u>=u0)[0] 								# Indices of electrons which have enough energy to cause excitation
		ll=np.where(np.random.rand(len(kk))<=p)[0]          # Generating a random vector and finding indices with value less than p where p is the probability of colllision
		kl=kk[ll]                   						# Indices of electrons which collided to give a photon  

		#P = np.random.rand(len(kl))						# generating random vector to reset position by after collision occurs 
		#xx[kl] = xx[kl]-dx[kl]*P             			 	# position reset by small factor with value less than the value of a single step
		#u[kl] = 0                                           # inelastic collision implies velocity reset to 0  

		
		"""
		This piece of code updates the values of electron position, velocity after collision considering the collision to happen at a time dt between (k-1)t and kt. A 
		random array of values of dt is generated. Position is changed by ((u[kl]-1)*dt + 0.5*(dt**2))+0.5*(1-dt)**2 and removing the default step dx. The acceleration 
		of electron is 1 unit and it accelerates for time (1-dt), hence, its velocity becomes 1-dt at the end of unit time step.
		"""

		dt = np.random.rand(len(kl))					
		xx[kl] = xx[kl]-dx[kl]+((u[kl]-1)*dt + 0.5*(dt**2))+0.5*(1-dt)**2
		u[kl] = 1-dt



		I.extend(xx[kl].tolist())							# Add the emitted phtons at those positions to intensity

		m = int(np.random.randn()*Msig + M)  				# Ejecting M number of electrons based around a normal distribution with mean M and stddev Msig						
		empty_xpos = np.where(xx==0)[0]						# Finding empty/available positions of electrons
		elec_generated = min(len(empty_xpos),m) 			# Adding the minimum of the empty positions and electrons ejected 
		xx[empty_xpos[0:elec_generated]] = 1.0				# Setting their positions as 1
		u[empty_xpos[0:elec_generated]] = 0.0				# Keeping velocity as 0 as they are ejected with 0 energy
		dx[empty_xpos[0:elec_generated]] = 0.0

		ii=np.where(xx>0)[0]								#Identifying present electrons
	
		X.extend(xx[ii].tolist())							#Extending position and velocity lists
		V.extend(u[ii].tolist()) 


def plot_electron_density():
	plt.figure(0)
	plt.hist(X,histtype='bar', bins=np.arange(0,n+1,n/100),ec='black',alpha=1) 
	plt.title('Electron Density')
	plt.xlabel(r'x$\rightarrow$',fontsize=12)
	plt.ylabel(r'Number of Electrons$\rightarrow$',fontsize=12)
	plt.show()

def plot_intensity():
	plt.figure(1)
	plt.hist(I,histtype='bar', bins=np.arange(0,n+1,n/100),ec='black',alpha=1) 
	plt.title('Emission intensity')
	plt.xlabel(r'x$\rightarrow$',fontsize=12)
	plt.ylabel(r'Intensity of Photons$\rightarrow$',fontsize=12)
	plt.show()

def plot_elec_phasespace():
	plt.figure(2)
	plt.plot(X,V,'x') 
	plt.title('Electron Phase Space')
	plt.xlabel(r'Position$\rightarrow$',fontsize=12)
	plt.ylabel(r'Velocity$\rightarrow$',fontsize=12)
	plt.show()

def tabulate_data():
	bins = plt.hist(I,bins=np.arange(0,n+1,n/100))[1]    	# Bin positions are obtained
	count = plt.hist(I,bins=np.arange(0,n+1,n/100))[0]   	# Population counts obtained
	xpos = 0.5*(bins[0:-1] + bins[1:]) 					#No. of bin end-points would be 1 greater, hence, the mean is used to get population of a particular bin
	df = pd.DataFrame()   								# A pandas dataframe is initialized to do the tabular plotting of values.
	df['Position'] = xpos
	df['count'] = count
	pd.set_option("display.max_rows", None, "display.max_columns", None)
	print(df)






#Driver Commands

#if arguments are given, checking if correct number of argumnets given
try:
	len(argv)==1 or len(argv)==7
except IOError:
	print("Invalid number of arguments. Usage example: python3 EE2703_ASSIGN6_EE19B077.py n=value M=value Msig=Value nk=Value u0=Value p=Value seed=Value ")
	exit(0)


#Update default values
if(len(argv)==7):
	n=int(argv[1])				
	M=int(argv[2])				
	nk=int(argv[3])
	u0=float(argv[4])
	p=float(argv[5])
	Msig=float(argv[6])    
	print("Using user provided parameters")

else:
	n=100                      	#spatial grid size                    
	M=5                        	#number of electrons injected per turn                    
	nk=500                     	#number of turns to simulate                    
	u0=5					   	#threshold velocity
	p=0.25					   	#probability of ionisation
	Msig=2                      #sigma of normal distribution

	print("Using defaults. If u want to use your own, put all 6 parameters as commandline arguments")

#initialising vectors for storing instantaneous information on all electrons
xx=np.zeros((n*M),dtype='float')	#electron postion
u=np.zeros((n*M),dtype='float')		#electron velocity
dx=np.zeros((n*M),dtype='float')	#displacement of current turn


#Initialising Lists
I=[]	#intensity of emitted light
X=[]	#position of electron
V=[]	#velocity of electron


electron_info()				#Simulates the tubelight using a loop

plot_electron_density()		#Plotting a population plot of electron density using hist
plot_intensity()			#Plotting a population plot of photon intensity using hist
plot_elec_phasespace()		#Plotting position of electrons vs their velocities
tabulate_data()				#Tabulating data for intensity vs position

