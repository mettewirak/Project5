import numpy as np
import matplotlib.pyplot as plt


def read(filename):
	with open(filename) as f:
		lines = f.readlines()
		energy = [line.split()[0] for line in lines]
		number_counted = [line.split()[1] for line in lines]
		probability = [line.split()[2] for line in lines]
	return energy, number_counted, probability
path="C:/Users/Mette Wirak/Documents/Faglig/Universitetet i Oslo/Computational Physics/build-Project5-Desktop_Qt_5_9_1_MinGW_32bit-Debug/" #Endre til din

filenames=["Histogram time=7, runs=3, N=500, dm=0.01.txt"]

fil_label=["$T=1$"]



i = 0

for fil in filenames:
	fil=path+fil
	energy, number_counted, probability=read(fil)
	plt.plot(energy, number_counted label=fil_label[i])
	i=+1


plt.xlabel('Amout of money a agent has [M0]')
plt.ylabel('Number of times the amount of money has appeared')

plt.legend(loc=2)
plt.show()
