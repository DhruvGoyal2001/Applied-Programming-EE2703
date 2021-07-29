"""
Assignment 7, EE2703 - Applied Lab Programming

Dhruv Goyal, EE19B077

Usage example: python3 EE2703_ASSIGN7_EE19B077.py               

The aim of this assignment is to analyze “Linear Time-invariant Systems” using the scipy.signal library in Python .We limit our analysis to systems with rational polynomial transfer functions.        
"""

# Importing libraries

import numpy as np 
import scipy.signal as sp 
import pylab as py
import matplotlib.pyplot as plt


#defining functions

def def_H(a,b):											# This function takes two values a,b and returns 2 1-D polynomials with shown coeffcients
	num = np.poly1d([1,a])								# Numerator
	den = np.poly1d([1,2*a, a**2 + b**2])				# Denominator
	return num,den

def build_cos(t,omega):									# This function takes time and frequency omega and returns an oscillating, decaying function of a spring system
	f1 = np.cos(omega*t)
	f2 = np.exp(-0.05*t)*np.heaviside(t,0.5)
	return f1*f2

def plot_response(t,x):									# Plots time response
	plt.plot(t,x, label = "Time response")
	plt.xlabel(r'time$\rightarrow$',fontsize=12)
	plt.ylabel(r'x(t)$\rightarrow$',fontsize=12)
	plt.legend(loc = "upper left")
	plt.grid(True)
	plt.show()

def plot_coupled(t1,t2,x,y):							# Plots time response x(t) and y(t) of a system of coupled equations		
	plt.figure(4)
	plt.plot(t1,x,"r",label = 'x(t)')
	plt.plot(t2,y,label = 'y(t)')
	plt.title(r"Solution for x(t) and y(t)")
	plt.xlabel(r't$\rightarrow$',fontsize=12)
	plt.ylabel(r'$x(t),y(t)\rightarrow$',fontsize=12)
	plt.legend(loc="upper right")
	plt.grid(True)
	plt.show()

def plot_RLC_bode(w,S,phi):								# Plots magnitude and phase bode plot of a given transfer function on semilog axis
	plt.figure(5)
	plt.subplot(2,1,1)
	plt.semilogx(w,S)
	plt.grid(True)
	plt.title("Magnitude and Phase plot of RLC circuit")
	plt.ylabel(r'$|H(s)|\rightarrow$')
	plt.xlabel(r'Frequency$\rightarrow$')
	plt.subplot(2,1,2)
	plt.semilogx(w,phi)
	plt.grid(True)
	plt.ylabel(r'$\angle(H(s))\rightarrow$')
	plt.xlabel(r'Frequency$\rightarrow$')
	plt.show()

def plot_rlc_response(i,y,time):						# Plots output time response when an input signal is given to a two port network with specific transfer function  
	plt.figure(i)
	plt.plot(T,y,label = "Time response for 0<time<{}".format(time))
	plt.xlabel(r't$\rightarrow$',fontsize=12)
	plt.ylabel(r'y(t)$\rightarrow$',fontsize=12)
	plt.legend(loc = "upper left")
	plt.grid(True)
	plt.show()


def build_v(T):											# Returns input signal vi(t), np.heaviside gives a unit step funciton with 0.5 being it's value at t=0, 1 for t>0 and 0 for t<0
	x = np.cos(1000*T)*np.heaviside(T,0.5)
	y = np.cos((10**6)*T)*np.heaviside(T,0.5)
	return x-y



# Driver Commands

num1,den1 = def_H(0.5,1.5)								# Returns num and den on the basis of given parameters
T = np.linspace(0,50,1000)								# Time array defined from o to 50s in 1000 steps 
den1 = np.polymul([1,0,2.25],den1) 						# Finding X(s)
H1 = sp.lti(num1,den1)									# Defining Transfer function
t,x = sp.impulse(H1,T=T)								# Using sp.impulse to find the impulse response of given transfer fucntion

plt.figure(1)
plt.title("Spring system with decay = 0.5")				# Plotting the time response
plot_response(t,x)

num2,den2 = def_H(0.05,1.5)								# Proceeding like before but with decay 0.05 instead of 0.5
den2 = np.polymul([1,0,2.25],den2) 
H2 = sp.lti(num2,den2)
t,x = sp.impulse(H2,T=T)

plt.figure(2)
plt.title("Spring system with decay = 0.05")
plot_response(t,x)


plt.figure(3)
i=0 
for alpha in np.arange(1.4,1.6,0.05):					# Changing frequency from 1.4 to 1.6 in steps of 0.05
	T = np.linspace(0,100,1000)							# Defining time array from 0 to 100s
	H = sp.lti([1],[1,0,2.25])							# Defining Transfer Funtion
	Cos = build_cos(T,alpha)							# Returns input signal
	t,y,svec = sp.lsim(H,Cos,T)							# Using sp.lsim to find convolution or the system response in time to the given input
	i+=1												# Plotting using subplots, labelpad defines distance between axis and it's label
	plt.subplot(3,2,i)									
	plt.xlabel(r't$\rightarrow$',fontsize=7,labelpad = -5)
	plt.ylabel(r'$Time Response\rightarrow$',fontsize=7,labelpad=-5)
	plt.plot(t,y)
plt.show()


T = np.linspace(0,20,1000)								# Defining time array from 0 to 20s
H_x = sp.lti(np.poly1d([1,0,2]),np.poly1d([1,0,3,0]))	# Defining the calculated transfer function of x and y using the given coupled equations and initial conditions of the system							
H_y = sp.lti(np.poly1d([2]),np.poly1d([1,0,3,0]))
t1,x = sp.impulse(H_x,T=T)								# Using sp.impulse to find the impulse response x(t), y(t)
t2,y = sp.impulse(H_y,T=T)
plot_coupled(t1,t2,x,y)									# Plotting the coupled time response 


H_rlc = sp.lti(np.poly1d([10**12]),np.poly1d([1,10**8,10**12])) # Transfer Function of the given two port network using given initial condition
w,S,phi=H_rlc.bode()											# Finding the bode plot magnitude and phase response using H.bode()
plot_RLC_bode(w,S,phi)											# Plotting the bode plots


T = np.linspace(0,30*(1e-6),1000)						# Defining Time array from 0 to 30us to analyze the effect of cos(10^6t) 
vi = build_v(T)											# Returns input signal array for given time array
t,y1,svec = sp.lsim(H_rlc,vi,T)							# Using sp.lsim to find the two port network response to given input
plot_rlc_response(6,y1,"30us")							# Plotting the time response


T = np.linspace(0,0.01,100000)							# Proceeding as before, just changing time definition upto 10ms to see the overall response
vi = build_v(T)											# The overall response is dominated by lower frequency component, cos(10^3t) as given circuit acts as a low pass filter
t,y2,svec = sp.lsim(H_rlc,vi,T)
plot_rlc_response(7,y2,"10ms")
