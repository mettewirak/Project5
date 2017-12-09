import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as Gamma
from pylab import *


def read(filename):
	with open(filename) as f:
		lines = f.readlines()
		energy = [line.split()[0] for line in lines]
		probability = [line.split()[2] for line in lines]		
	return energy, probability

path="C:/Users/Mette Wirak/Documents/Faglig/Universitetet i Oslo/Computational Physics/5 - Histogram data/"

filenames=["Histogram N=500 time=7 runs=4 lambda=0.000000 alpha=0.000000.txt", "Histogram N=500 time=7 runs=4 lambda=0.250000 alpha=0.000000.txt", "Histogram N=500 time=7 runs=4 lambda=0.500000 alpha=0.000000.txt"]
fil_label=["lambda = 0.00", "lambda = 0.25", "lambda = 0.50", "lambda = 0.90"]

i = 0
for fil in filenames:
	fil=path+fil
	energy, probability=read(fil)

	plt.plot(energy, probability, label=fil_label[i])
	i=i+1


# Finding the analytic
m = np.arange(0., 10., 0.01)
plt.plot(m, (np.power(1,1)/Gamma(1))*np.power(m,0)*np.exp(-1*m), '--')
plt.plot(m, (np.power(2,2)/Gamma(2))*np.power(m,1)*np.exp(-2*m), '--')
plt.plot(m, (np.power(4,4)/Gamma(4))*np.power(m,3)*np.exp(-4*m), '--')
#plt.plot(m, (np.power(28,28)/Gamma(28))*np.power(m,27)*np.exp(-28*m), '--', label="analytic lambda = 0.90")


plt.xlabel('log(m)')
plt.ylabel('log(f(m))')

plt.xscale('log')
plt.yscale('log')

plt.legend(loc=1)
plt.show()
