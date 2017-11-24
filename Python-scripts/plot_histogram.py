import numpy as np
import matplotlib.pyplot as plt


def read(filename):
	with open(filename) as f:
		lines = f.readlines()
		energy=np.zeros(len(lines))
		number_counted=np.zeros(len(lines))
		probability=np.zeros(len(lines))
		i=0
		for line in lines:
			words= line.split()

			energy[i]=words[0];
			number_counted[i]=words[1]
			probability[i]=words[2]
			i+=1

	return energy, number_counted,probability
#path="/home/arnlaug/Documents/UiO/Compfys/Project/5/" #Endre til din
path="/home/arnlaug/Documents/UiO/Compfys/build-5_Arnlaug-Desktop-Release/"


#filenames=["Histogram time=7, runs=1, N=500, dm=0.01.txt"]
filenames=["Test.txt", "Test2.txt"]
#"Histogram time=7, runs=1, N=500, dm=0.01.txt"
fil_label=["$runs=10^0$","$runs=10^1$" ]



i = 0

for fil in filenames:
	fil=path+fil
	energy, number_counted, probability=read(fil)

	plt.plot(energy, probability, label=fil_label[i])
	i=+1


plt.xlabel('Amout of money a agent has [M0]')
plt.ylabel('Number of times the amount of money has appeared')

plt.legend(loc=2)
plt.show()
