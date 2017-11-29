import numpy as np
import matplotlib.pyplot as plt


def read(filename):
	with open(filename) as f:
		lines = f.readlines()
		energy = [line.split()[0] for line in lines]
		probability = [line.split()[2] for line in lines]		
	return energy, probability

path="/uio/hume/student-u87/mettewir/Documents/Computational Physics/Project5/Histogram data/"

filenames=["Histogram N=500 time=7 runs=4 lambda=0.000000 alpha=2.000000.txt", "Histogram N=500 time=7 runs=4 lambda=0.000000 alpha=1.500000.txt", "Histogram N=500 time=7 runs=4 lambda=0.000000 alpha=1.000000.txt", "Histogram N=500 time=7 runs=4 lambda=0.000000 alpha=0.500000.txt", "Histogram N=500 time=7 runs=4 lambda=0.000000 alpha=0.000000.txt", "Histogram N=1000 time=7 runs=4 lambda=0.000000 alpha=2.000000.txt", "Histogram N=1000 time=7 runs=4 lambda=0.000000 alpha=1.500000.txt", "Histogram N=1000 time=7 runs=4 lambda=0.000000 alpha=1.000000.txt", "Histogram N=1000 time=7 runs=4 lambda=0.000000 alpha=0.500000.txt", "Histogram N=1000 time=7 runs=4 lambda=0.000000 alpha=0.000000.txt"]

fil_label=["l = 0,a = 2.0","l = 0, a = 1.5", "l = 0, a = 1.0", "l = 0, a = 0.5", "l = 0, a = 0.0", "l = 0, a = 2.0","l = 0, a = 1.5", "l = 0, a = 1.0", "l = 0, a = 0.5", "l = 0, a = 0.0"]


i = 0
for fil in filenames:
	fil=path+fil
	energy, probability=read(fil)

	plt.plot(energy, probability, label=fil_label[i])
	i=i+1


plt.xlabel('m')
plt.ylabel('f(m)')

plt.legend(loc=1)
plt.show()
