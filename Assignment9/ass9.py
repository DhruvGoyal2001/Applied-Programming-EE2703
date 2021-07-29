"""
Assignment 9, EE2703 - Applied Lab Programming

Dhruv Goyal, EE19B077

Usage example: python3 EE2703_ASSIGN9_E19B077.py               

In this assignment, we explore how to obtain the DFT, and how to recover theanalog Fourier Tranform for some known functions by the proper sampling of the function.      
"""

# Importing libraries
from pylab import *

#Defining dictionaries for titles in graphs
title_dict = {'sin': r"Spectrum of $sin(5*t)$",'cos': r"Spectrum of $1 + 0.1cos(t))cos(10t)$", 'cos3': r"Spectrum of $sin^3(t)$",'sin3': r"Spectrum of $cos^3(t)$",'fm': r"Spectrum of $cos(20t+5cos(t))$", 'gauss': r"Spectrum of $\exp(-t^2/2)$"}

# Function to find fft
def fft_(f, N, func, xLim = 15):
	t = linspace(-N*pi/128, N*pi/128, N + 1)[:-1]			# Defining time as an array
	y = f(t)												# Defining y as the function of t			
	Y = fftshift(fft(y))/ float(N)							# Using fftshift on the fft of function y
	w = linspace(-64, 64, N + 1)[:-1]						# Defining frequency as an array
	figure()												# Plotting Magnitude and Phase spectrum
	subplot(2, 1, 1)
	plot(w, abs(Y), lw = 2, label = "Magnitude")
	ttl = title_dict[func]
	title(ttl)
	ylabel(r"$|Y|\rightarrow$", fontsize = 12)
	legend(loc = "best")
	xlim([-xLim, xLim])										# Limiting values on x-axis
	grid(True)
	subplot(2, 1, 2)
	if func != 'fm':										
		plot(w, angle(Y),'ro', markersize = 3, label = "Phase")
	ylabel(r"Phase of $Y\rightarrow$", fontsize=12)
	xlabel(r"$\omega\rightarrow$", fontsize=12)
	ii = where(abs(Y) > 1e-3)
	plot(w[ii], angle(Y[ii]), 'go', label = "Phase at peaks")
	legend(loc = "upper left")
	if func == 'fm':
		legend(loc = "center")
	xlim([-xLim, xLim])
	grid(True)
	show()


def estimate(N, T):											# Defining function to calculate spectrum of Gaussian function 
	t = linspace(- T/2, T/2, N + 1)[:-1]
	w = linspace(- N*pi/T, N*pi/T, N + 1)[:-1]
	y = exp(-0.5 * t**2)
	Y_true = exp(-0.5 * w**2)*sqrt(2*pi)
	Y = fftshift(abs(fft(y)))/N
	Y = Y*sqrt(2*pi)/max(Y)

	return max(abs(Y - Y_true)), w, Y, Y_true



# Driver Commands

# FFT for different functions
fft_(lambda t : sin(5*t), 128, 'sin')
fft_(lambda t : (1 + 0.1*cos(t))*cos(10*t), 512, 'cos')
fft_(lambda t : sin(t)**3, 512, 'cos3')
fft_(lambda t : cos(t)**3, 512, 'sin3')
fft_(lambda t : cos(20*t + 5*cos(t)), 512, 'fm', xLim = 35)


# Using while loop to iterate 
i = 1
while estimate(N = 512, T = i*pi)[0] > 1e-15:
	i += 1

print('Time range for accurate spectrum : ' + str(i) + 'pi')
print('Accuracy : ' + str(estimate(N = 512, T = i * pi)[0]))

w, Y, Y_true = estimate(N = 512, T = i*pi)[1:]


# Plotting Gaussian function Magnitude and Phase spectrum 
xLim = 10
print(len(Y))
figure()
subplot(2, 1, 1)
plot(w, abs(Y), lw = 2, label = "Magnitude")
legend(loc = "best")
ylabel(r"$|Y|\rightarrow$", fontsize = 12)
xlim([-xLim, xLim])
grid(True)
subplot(2, 1, 2)
plot(w, angle(Y),'ro', markersize = 3, label = "Phase")
ii = where(abs(Y) > 1e-3)
plot(w[ii], angle(Y[ii]), 'go', label = "Phase at peaks")
ylabel(r"Phase of $Y\rightarrow$", fontsize=12)
xlabel(r"$\omega\rightarrow$", fontsize=12)
xlim([-xLim, xLim])
legend(loc = "best")
show()

