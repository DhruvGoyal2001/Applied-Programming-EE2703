"""
Assignment 10, EE2703 - Applied Lab Programming

Dhruv Goyal, EE19B077

Usage example: python3 EE2703_ASSIGN10_EE19B077.py               

In this assignment, we explore how to obtain the DFT and analyze non-periodic functions as well      
"""

# Importing libraries

from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm
import numpy as np


# Defining functions

# Defining dictionaries for title and functions
title_dict = {'sin': r"Spectrum of $\sin(\sqrt{2}t)$", 'cos3' : r"Spectrum of $\cos^3(0.86t)$",'cos':r"Spectrum of $\cos(\omega*t + \delta)$",'chirp':r"Spectrum of $\cos\left(16\left(1.5+\frac{t}{2\pi}\right)t\right)$"}
func_dict = {'sin':lambda x: sin(sqrt(2)*x),'cos3':lambda x: cos(0.86*x)**3,'cos':lambda x: cos(omega*x + delta),'chirp' : lambda t : cos(16*t*(1.5+t/(2*pi)))}


#Function to find a plot the FFT spectrum of various functions
def find_fft(func,N,t_l=-pi,t_r=pi,windowing = False,xLim=8, noise = False,Plot=True):

    t = linspace(t_l,t_r,N+1); t=t[:-1]                         # Defining time array
    
    fmax = 1/(t[1]-t[0])                                        # Defining fmax and importing function from dictionary
    y = func_dict[func](t)
    
    if windowing == True:                                       # Condition for multiplying by hamming window or not   
        y=y*hamming(N)

    if noise == True:                                           # Condition for adding noise
        y = y + 0.1*randn(N)

    y[0]=0                                                      # Calculating the FFT and using fftshift for symmetric
    y=fftshift(y)
    Y=fftshift(fft(y))/float(N)
    w=linspace(-pi*fmax,pi*fmax,N+1);w=w[:-1]

    if Plot == True:
        figure()                                                # Plotting and labelling Magnitude spectrum
        subplot(2, 1, 1)
        plot(w, abs(Y), lw = 2, label = "Magnitude")
        ttl = title_dict[func]                                  # Extracting title from the dictionary
        title(ttl)
        ylabel(r"$|Y|\rightarrow$", fontsize = 12)
        legend(loc = "best")
        xlim([-xLim, xLim])                                     # Limiting values on x-axis
        grid(True)
    
        subplot(2, 1, 2)                                        # Plotting and labelling phase spectrum
        plot(w, angle(Y),'ro', markersize = 4, label = "Phase")
        ylabel(r"Phase of $Y\rightarrow$", fontsize=12)
        xlabel(r"$\omega\rightarrow$", fontsize=12)
        legend(loc = "best")
        xlim([-xLim, xLim])
        grid(True)
        show()

    if Plot == False:                                           # Returning Y and w if not to be plotted
        return Y,w


#Function to return value of Hamming window for a given length N
def hamming(N):
    n = arange(N)
    return fftshift(0.54+0.46*(np.cos(2*pi*n/(N-1))))


#Function to plot graph of Sin(1.414t) in a periodic manner
def plot_sin():
    figure(1)                                                   # Without Hamming window
    t=np.linspace(-pi,pi,64)
    plot(t,y_sin(t),'ro', label = "") 
    plot(t+2*pi,y_sin(t),'bo') 
    plot(t-2*pi,y_sin(t),'bo')
    ylabel(r"$\sin(\sqrt{2}t)\rightarrow$", fontsize=12)
    xlabel(r"$Time\rightarrow$", fontsize=12)
    title('Without hamming window')
    grid(True)
    show()


    figure(2)                                                    # With Hamming window
    plot(t,y_sin_ham(t),'ro')
    plot(t+2*pi,y_sin_ham(t),'bo')
    plot(t-2*pi,y_sin_ham(t),'bo')
    ylabel(r"$\sin(\sqrt{2}t)*(Hamming Window)\rightarrow$", fontsize=12)
    xlabel(r"$Time\rightarrow$", fontsize=12)
    title('With hamming window')
    grid(True)
    show() 


#Function to estimate the values of fequency and phase for a given unknown cosine function
def estimate(w,y):
    ind = where(w>0)[0]                                         # Dividing into two halves - +ve and -ve
    y = y*hamming(128)                                          # Finding FFT of y, windowed
    y[0]=0
    y = fftshift(y)
    Y = fftshift(fft(y))/128.0                                  
    w0 = sum((abs(Y[ind])**2)*w[ind])/sum(abs(Y[ind])**2)       # Using weighted average to estimate frequency
    delta = angle(Y[::-1][argmax(abs(Y[::-1]))])                # Estimating delta by finding the angle at the peaks or maximum of |Y|
    return w0,delta


#Function to define and return the time and frequency array
def t_w():
    t = (linspace(-2*pi,2*pi,129))[:-1]
    dt = t[1]-t[0]
    fmax = 1/dt
    w = (linspace(-pi*fmax,pi*fmax,129))[:-1]
    return t,w


#Function to plot the chirped signal as a function of time
def plot_chirp():
    time = linspace(-pi,pi,1024)
    plot(time,func_dict['chirp'](time))
    ylabel(r"$\cos\left(16\left(1.5+\frac{t}{2\pi}\right)t\right)\rightarrow$", fontsize=12)
    xlabel(r"$Time\rightarrow$", fontsize=12)
    title('Chirp Signal')
    grid(True)
    show()


#Function to calculate the FFT of chirped signal after breaking the array into 16 parts
def chirp_broken():
    Y = []
    for i in range(16):
        t_left = -pi + (2*pi)*i/16
        t_right = -pi + (2*pi)*(i+1)/16
        Y.append(find_fft('chirp',64,t_left,t_right,windowing=False,Plot=False)[0])
    return np.array(Y)


#Function to plot the surafce plots for the calculated FFT for each of the parts of chirp signal
def plot_surface():
    Y_arr = chirp_broken()
    Y_arr = Y_arr.T

    t1 = linspace(-pi,pi,16)
    t = (linspace(-pi,pi,1025))[:-1]
    fmax=1/(t[1]-t[0])
    w = (linspace(-pi*fmax,pi*fmax,65))[:-1]

    t1,w = meshgrid(t1,w)
    print(t1.shape)

    ax = p3.Axes3D(figure())
    surf = ax.plot_surface(w,t1,abs(Y_arr),rstride=1,cstride=1,cmap=cm.jet)
    colorbar(surf, shrink = 0.4, aspect = 5)
    ylabel(r"$\omega\rightarrow$", fontsize=12)
    xlabel(r"$Time\rightarrow$", fontsize=12)
    title('Surface time-frequency plot')
    ax.set_zlabel(r'$|Y|\rightarrow$')
    show()



# Driver commands

y_sin=lambda t: sin(sqrt(2)*t)                                          # Defining sin and windowed sin for plotting v/s time                           
y_sin_ham = lambda t: (sin(sqrt(2)*t))*hamming(len(t))                  

plot_sin()                                                              # Plotting sin(sqrt(2)t) v/s time in periods of 2*pi

find_fft('sin',256,-4*pi,4*pi,windowing=False)                          # Calculating and plotting FFT for sin(sqrt(2)t)
find_fft('sin',256,-4*pi,4*pi,True)                                     # Repeating the same but using hamming window on the function 

find_fft('cos3',256,-4*pi,4*pi,False)                                   # Repeating the process for cos^3(0.86t)
find_fft('cos3',256,-4*pi,4*pi,True)


omega, delta = np.random.uniform(0.5,1.5), np.random.uniform(0, pi)     # Using random.uniform to get input values of omega and delta 
t,w = t_w()                                                             # Defines time and frequency arrays

find_fft('cos',128,-pi,pi,False,8,False)                                # Spectrum for cos(omega*t + delta) is plotted using above values
W,Y = estimate(w, cos(omega*t + delta))                                 # Estimate function is used to determine the frequency and phase
print("Input frequency = {}, delta = {}".format(omega,delta))
print("Estimated frequency = {}, delta = {}".format(W,Y))

find_fft('cos',128,-pi,pi,False,8,True)                                 # Repeating for signal with added white Gaussian noise
W1,Y1 = estimate(w, cos(omega*t + delta) + 0.1*randn(128))              # Noise generated using randn function
print("With noise, Estimated frequency = {}, delta = {}".format(W1,Y1))

plot_chirp()                                                            # Chirped signal is plotted against time
find_fft('chirp',1024,-pi,pi,windowing=False,xLim=50)                   # Spectrum for chirped signal is evaluated and plotted

plot_surface()                                                          # Surface plots for broken "chirped" signal are obtained

