import numpy as np
import matplotlib.pyplot as plt


def read(filename):
	with open(filename) as f:
		lines = f.readlines()
		energy = [line.split()[0] for line in lines]
		probability = [line.split()[2] for line in lines]		
	return energy, probability

path="C:/Users/Mette Wirak/Documents/Faglig/Universitetet i Oslo/Computational Physics/5 - Histogram data/"

filenames=["Histogram N=500 time=7 runs=4 lambda=0.500000 alpha=0.500000.txt", "Histogram N=500 time=7 runs=4 lambda=0.500000 alpha=1.000000.txt", "Histogram N=500 time=7 runs=4 lambda=0.500000 alpha=1.500000.txt", "Histogram N=500 time=7 runs=4 lambda=0.500000 alpha=2.000000.txt"]

fil_label=["alpha = 0.5", "alpha = 1.0", "alpha = 1.5", "alpha = 2.0"]

i = 0
for fil in filenames:
	fil=path+fil
	energy, probability=read(fil)

	plt.plot(energy, probability, label=fil_label[i])
	i=i+1


#m = np.arange(0., 10., 0.01)
#plt.plot(m, np.exp(-m), '--', label="$Analytical$")


plt.xlabel('m')
plt.ylabel('f(m)')

plt.xscale('log')
plt.yscale('log')

plt.legend(loc=1)
plt.show()
