"""
Assignment 8, EE2703 - Applied Lab Programming

Dhruv Goyal, EE19B077

Usage example: python3 EE2703_ASSIGN8_EE19B077.py               

The aim of this assignment is to analyze circuits using symp library in Python .We limit our analysis to systems with rational polynomial transfer functions.        
"""

# Importing libraries
import numpy as np
import sympy as sym
import scipy.signal as sp
import pylab as py
import matplotlib.pyplot as plt

# Defining Functions


def lowpass(R1,R2,C1,C2,G,Vi): 	# Lowpass filter takes circuit parameters
    s=sym.symbols('s')																				# Defining s as symbol
    A=sym.Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])	# Nodal analysis parameter matrix
    b=sym.Matrix([0,0,0,-Vi/R1])																	# independent source matrix
    V=A.inv()*b																						# Obtaning node voltages
    return A,b,V


def mag_response(Vo):
	s=sym.symbols('s')
	w = np.logspace(0,8,801)																		# Defining frequency array
	ss=1j*w																							# Converting to jw form																						
	hf=sym.lambdify(s,Vo,'numpy')																	# returns temporary lamda function in python from sympy funciton for vector evaluation
	v=hf(ss)
	return v,w 																						# returning magnitude and frequency


def convolve(Vo,T,f):																				# Function to perform convolution using sp.lsim
	num, den = symToH(Vo)																			# symToH breaks into numerator and denominator
	H = sp.lti(num,den)
	t,y,svec = sp.lsim(H,f,T)
	return t,y


def highpass(R1,R2,C1,C2,G,Vi):		# Highpass Filter
    s=sym.symbols('s')
    A=sym.Matrix([[0,0,1,-1/G],[-1/(1+1/(s*R2*C2)),1,0,0],[0,-G,G,1],[-s*C1-s*C2-1/R1,s*C2,0,1/R1]])
    b=sym.Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return A,b,V

def symToH(H):
    #H = sym.simplify(H)																			# Simplifies the given expression
    n,d = sym.simplify(H).as_numer_denom()															# Returns numerator and denominator
    n,d = sym.poly(n,s), sym.poly(d,s)																# Transforms expression into a polynomial
    num,den = n.all_coeffs(), d.all_coeffs()														# Returns coeffcients
    num,den = [float(f) for f in num], [float(f) for f in den]										# Converting all to floats
    return num,den 																					# Returns numerator and denominator


def damped_sin(t, freq, decay):																		# Returns damped sinusoid funtion
	return np.cos(2*np.pi*freq*t)*np.exp(-decay*t)

# Driver Commands

s=sym.symbols('s')

# Finding Lowpass filter impulse response
A,b,V_il=lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vo_il=V_il[3]												# Laplace transform of Vo for impulse
v_il,w_il=mag_response(Vo_il)								# Magnitude of impulse response


# Finding Lowpass filter step response
A,b,V_ul=lowpass(10000,10000,1e-9,1e-9,1.586,1/s)			# Solving for input as unit step signal
Vo_ul=V_ul[3]												# Laplace transform of Vo for step input
v_ul,w_ul=mag_response(Vo_ul)								# Magnitude of step response


# Finding Highpass filter impulse response
A,b,V_ih=highpass(10000,10000,1e-9,1e-9,1.586,1)
Vo_ih=V_ih[3]												# Laplace transform of Vo for impulse
v_ih,w_ih=mag_response(Vo_ih)								# Magnitude of impulse response


# Finding Highpass filter step response
A,b,V_uh=highpass(10000,10000,1e-9,1e-9,1.586,1/s)			# Solving for input as unit step signal
Vo_uh=V_uh[3]												# Laplace transform of Vo for step input
v_uh,w_uh=mag_response(Vo_uh)								# Magnitude of step response


# Finding step response
T = np.linspace(0,0.005,10000)								# Defining time array
u=np.heaviside(T,1)											# defines step function
t_ul, y_ul = convolve(Vo_il,T,u)							# Finding step response by defining the impulse response and convolving for lowpass
#t_uh, y_uh = convolve(Vo_ih,T,u)							# Finding step response for high pass filter


# Finding response of the below defined input
T=np.linspace(0,0.001,100000)
Sin=(np.sin(2000*np.pi*T)+np.cos(2*np.pi*T*1e6))

t_Sinl,y_Sinl=convolve(Vo_il,T,Sin)							# For lowpass filter
t_Sinh,y_Sinh=convolve(Vo_ih,T,Sin)							# For highpass filter

u=np.heaviside(T,1)											# defines step function
t_uh, y_uh = convolve(Vo_ih,T,u)							# Finding step response for high pass filter

# Finding response to decaying sinusoid with high frequency
dSin=damped_sin(T,1e7,3000)									# Decaying sinusoid function

t_dsinl,y_dsinl=convolve(Vo_il,T,dSin)						# For lowpass filter
t_dsinh,y_dsinh=convolve(Vo_ih,T,dSin)						# For Highpass filter



# Finding response to decaying sinusoid with low frequency
T1 = np.linspace(0,0.5,1000)
dSin2=damped_sin(T1,1e2,50)									# Decaying sinusoid function

t_dsinl2,y_dsinl2=convolve(Vo_il,T1,dSin2)						# For lowpass filter
t_dsinh2,y_dsinh2=convolve(Vo_ih,T1,dSin2)						# For Highpass filter


#Plotting all the graphs

plt.figure(1)
plt.title("Magnitude response for lowpass filter")
plt.loglog(w_il, abs(v_il), lw=2, label = "Magnitude response for lowpass filter")
plt.xlabel(r'Frequency$\rightarrow$',fontsize=12)
plt.ylabel(r'Magnitude$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(2)
plt.title("Magnitude response for highpass filter")
plt.loglog(w_ih, abs(v_ih), lw=2, label = "Magnitude response for highpass filter")
plt.xlabel(r'Frequency$\rightarrow$',fontsize=12)
plt.ylabel(r'Magnitude$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(3)
plt.title("Unit step response for lowpass filter")
plt.plot(t_ul, y_ul, label = "Unit step response for lowpass filter")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vout$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(4)
plt.title("Unit step time response for highpass filter")
plt.plot(t_uh, y_uh, label = "Unit step response")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vout$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(5)
plt.title("Sum of sinusoids")
plt.plot(T, Sin, label = "Input voltage signal")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vout$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()


plt.figure(6)
plt.title("Time response for sum of sinusoids of highpass filter")
plt.plot(t_Sinh, y_Sinh, label = "Output Voltage for sum of sinusoids ")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vout$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(7)
plt.title("Time response for sum of sinusoids of lowpass filter")
plt.plot(t_Sinl, y_Sinl, label = "Output Voltage for sum of sinusoids ")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vout$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(8)
plt.title("High frequency damped sinusoid")
plt.plot(T, dSin, label = "Vin")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vin$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(9)
plt.title("Time response of decaying sinusoid with high frequency of lowpass filter")
plt.plot(t_dsinl, y_dsinl, label = "Output Voltage for decaying sinusoid ")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vout$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(10)
plt.title("Time response of decaying sinusoid with high frequency of highpass filter")
plt.plot(t_dsinh, y_dsinh, label = "Output Voltage for decaying sinusoid ")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vout$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(11)
plt.title("Low frequency damped sinusoid")
plt.plot(T1, dSin2, label = "Vin")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vin$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(12)
plt.title("Time response of decaying sinusoid with low frequency of lowpass filter")
plt.plot(t_dsinl2, y_dsinl2, label = "Output Voltage for decaying sinusoid ")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vout$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()

plt.figure(13)
plt.title("Time response of decaying sinusoid with low frequency of highpass filter")
plt.plot(t_dsinh2, y_dsinh2, label = "Output Voltage for decaying sinusoid ")
plt.xlabel(r'Time$\rightarrow$',fontsize=12)
plt.ylabel(r'Vout$\rightarrow$',fontsize=12)
plt.legend(loc = "best")
plt.grid(True)
plt.show()