"""
Assignment 3, EE2703 - Applied Lab Programming

Dhruv Goyal, EE19B077

Usage example: python3 EE2703_ASSIGN3_EE19B077.py
"""

# importing libraries
import matplotlib.pyplot as plt
import scipy
import scipy.special as sp
import pylab
import numpy as np

# noise variations
sigma=np.logspace(-1,-3,9)

#plotting the different noise level graphs and true graph
def noise_with_true_value(x,y):
    plt.figure(0,figsize=(7,6))
    plt.title(r'Plot with Added Noise levels')
    for i in range(9):                                          #all sets of values of y
        plt.plot(x,y[:,i],label=r'$\sigma$=%.3f'%sigma[i])
    plt.plot(x,g(x,1.05,-0.105),'-k',label='True curve')        #plotting the true curve given by exact bessel function
    plt.legend(loc='upper right')                               
    plt.ylabel(r'Function + Noise',fontsize = 14)
    plt.xlabel(r'time',fontsize = 14)
    plt.grid(True)
    plt.show()
    return

def plot_with_noise(x,y):
    plt.figure(0,figsize=(7,6))
    plt.title(r'Plot with Added Noise levels')
    for i in range(9):                                          #all sets of values of y
        plt.plot(x,y[:,i],label=r'$\sigma$=%.3f'%sigma[i])
    plt.legend(loc='upper right')                               
    plt.ylabel(r'Function + Noise',fontsize = 14)
    plt.xlabel(r'time',fontsize = 14)
    plt.grid(True)
    plt.show()
    return

#defining bessel function for time instant t
def g(t,A=1.05,B=-0.105):
    return A*sp.jn(2,t)+B*t

#Plotting the error bar 
def plot_error(x,y,i):
    y_true = g(x,1.0,-0.105)
    dif = np.std(y[:,i]-y_true)
    plt.plot(x,y_true,'b',label='true')                                           #plotting original curve
    plt.title('Q.5. Plot for sigma = ' + str(sigma[i]) + ' with error bars')
    plt.xlabel(r't$\rightarrow$',fontsize=14)
    plt.errorbar(x[::5],y[::5,i],dif,fmt='ro',label='errorbar')                       #plotting Error bars for every fifth point
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.show()
    return

# matrix which has parameters of the function
def generateM(x):
    fn_column = sp.jn(2,x)
    M = np.c_[fn_column,x]
    return M


# Calculating mean squared error
def find_error_matrix(x,y):
    E = np.zeros((21,21,9))
    A = np.linspace(0,2,21)
    B = np.linspace(-0.2,0,21)

    for k in range(9):
        f = y[:,k]
        for i in range(21):
            for j in range(21):
                E[i][j][k] = np.square((f -np.array(g(x,A[i],B[j])))).mean()

    return E

# plotting the mean squared error on a contour plot
def plot_contour(x,y):
    A = np.linspace(0,2,21)
    B = np.linspace(-0.2,0,21)
    error = find_error_matrix(x,y)
    plot = plt.contour(A,B,error[:,:,0],20)   
    plt.clabel(plot,inline=1,fontsize=10)
    plt.title('Contour Plot of error')
    plt.xlabel(r'$A$',size=10)
    plt.ylabel(r'$B$',size=10)

    # Finding minima of the function
    a = np.unravel_index(np.argmin(error[:,:,0]),error[:,:,0].shape)
    plt.plot(A[a[0]],B[a[1]],'o',markersize=3)
    plt.annotate('(%0.2f,%0.2f)'%(A[a[0]],B[a[1]]),(A[a[0]],B[a[1]]))   # Labeling the point of minimum constrained by values of A and B
    plt.show()
    return


# Testing equality by Matrix multiplication 
def testing_if_equal(x,C):
    X = generateM(x)
    H1 = np.array(g(x))
    H = np.matmul(X,C)
    return np.array_equal(H,H1)           #return true for equal

#Estimating A and B 
def estimateAB(M,b):
    return scipy.linalg.lstsq(M,b)

# Plotting Error vs standard deviation
def plot_errorAB(x,y):
    error_1 = np.zeros(9)
    error_2 = np.zeros(9)
    for i in range(9):                                         #predicting A and B and the error associated
        prediction,_,_,_ = estimateAB(generateM(x),y[:,i])
        error_1[i],error_2[i] = np.abs(prediction[0]-1.05),np.abs(prediction[1]-(-0.105))
    # Plotting error
    plt.plot(sigma,error_1,'r--',label='Aerr')
    plt.scatter(sigma,error_1)
    plt.plot(sigma,error_2, 'b--',label='Berr')
    plt.scatter(sigma,error_2)
    plt.legend(loc='upper left')
    plt.title("Variation Of error with Noise")
    plt.ylabel(r'MS Error$\rightarrow$',fontsize=15)
    plt.xlabel(r'$\sigma_{n}\rightarrow$',fontsize=15)
    plt.grid(True)
    plt.show()

    #plotting error in A and B on log Scale
    plt.loglog(sigma,error_1,'ro',label='Aerr')
    plt.stem(sigma,error_1,'-ro')
    plt.loglog(sigma,error_2,'bo',label='Berr')
    plt.stem(sigma,(error_2),'-bo')
    plt.legend(loc='upper left')
    plt.title("Variation Of error with Noise on loglog scale")
    plt.xlabel(r'$\sigma_{n}\rightarrow$',fontsize=15)
    plt.ylabel(r'MS Error$\rightarrow$',fontsize=15)
    plt.grid(True)
    plt.show()

    return


#loading the data from the generated file
data = np.loadtxt("fitting.dat")
time = data[:,0]                     #time variable
y = data[:,1:]                    #value of y for different noises

plot_with_noise(time,y)

noise_with_true_value(time,y)

plot_error(time,y,0)

# Creating M and parameter array, checking for equality
A = 1.05; B = -0.105
C = np.array([A,B]) 
print("Equality of M = ",testing_if_equal(time,C))

plot_contour(time,y)

prediction,_,_,_ = estimateAB(generateM(time),y[:,1])

print("Prediction of A,B= ",prediction)
print("Mean Squared Error in the A,B = ",np.abs(prediction[0]-1.05),np.abs(prediction[1]-(-0.105)))

#Plotting error vs stdev on linear and log scale
plot_errorAB(time,y)

