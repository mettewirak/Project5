import numpy as np
import matplotlib.pyplot as plt


def read(filename):
	with open(filename) as f:
		lines = f.readlines()
		amount=np.zeros(len(lines))
		number_counted=np.zeros(len(lines))
		probability=np.zeros(len(lines))
		log=np.zeros(len(lines))
		i=0
		for line in lines:
			words= line.split()

			amount[i]=words[0];
			number_counted[i]=words[1]
			probability[i]=words[2]
			log[i]=words[3]
			i+=1


	return amount, number_counted,probability,log


def read_gjennomsnitt(filename):
	with open(filename) as f:
		lines = f.readlines()
		agent=np.zeros(len(lines))
		amount=np.zeros(len(lines))

		i=0
		for line in lines:

			words= line.split()

			agent[i]=words[0];
			amount[i]=words[1]

			i+=1


	return agent,amount
#path="/home/arnlaug/Documents/UiO/Compfys/Project/5/" #Endre til din
path="/home/arnlaug/Documents/UiO/Compfys/build-5_Arnlaug-Desktop-Release/"


#filenames=["Gjennomsnittsformuer N=1000 time=7 runs=4 gamma=0 lambda=0 alpha =1.txt","Gjennomsnittsformuer N=1000 time=7 runs=4 gamma=1 alpha=2 lambda=0.txt","Gjennomsnittsformuer N=1000 time=7 runs=4 gamma=2 alpha=1 lambda=0.txt","Gjennomsnittsformuer N=1000 time=7 runs=4 gamma=3 alpha=1 lambda=0.txt","Gjennomsnittsformuer N=1000 time=7 runs=4 gamma=4 alpha=1 lambda=0.txt"]
#filenames=["Histogram N=1000 time=7 runs=4 gamma=0 lambda=0 alpha =1.txt","Histogram N=1000 time=7 runs=4 gamma=1 alpha=2 lambda=0.txt","Histogram N=1000 time=7 runs=4 gamma=2 alpha=1 lambda=0.txt","Histogram N=1000 time=7 runs=4 gamma=3 alpha=1 lambda=0.txt","Histogram N=1000 time=7 runs=4 gamma=4 alpha=1 lambda=0.txt"]
filenames=["Histogram N=1000 time=7 runs=4 gamma=1 alpha=1 lambda=0.txt","Histogram N=1000 time=7 runs=4 gamma=4 alpha=1 lambda=0.txt"]#,"Histogram N=1000 time=7 runs=4 gamma=3 alpha=1 lambda=0.txt","Histogram N=1000 time=7 runs=4 gamma=4 alpha=1 lambda=0.txt"]


#filenames=["Histogram N=500 time=7 runs=1 gamma=1 alpha=2 lambda=0.txt","Histogram N=500 time=7 runs=1 gamma=2 alpha=2 lambda=0.txt","Histogram N=500 time=7 runs=1 gamma=3 alpha=2 lambda=0.txt","Histogram N=500 time=7 runs=1 gamma=4 alpha=2 lambda=0.txt"]
fil_label=["0","1","2","3","4"]
color=['r','b']
i = 0

for fil in filenames:
	print fil
	fil=path+fil
	amount, number_counted, probability, log=read(fil)
	#agent,amount=read_gjennomsnitt(fil)
	fil_label[i]=plt.bar(amount, number_counted, label=fil_label[i], width=0.005,edgecolor = "none", color=color[i])
	#plt.plot( amount, number_counted, label=fil_label[i])
	#plt.plot(amount,np.exp(-amount), label=fil_label[i])
	#plt.plot(np.log10(agent), np.log10(amount), label=fil_label[i])
	i+=1
	



plt.xlabel('Amout of money a agent has [M0]')
plt.ylabel('Number of times the amount of money has appeared')

plt.legend(loc=2)
#plt.savefig('Sammenligne for runs')
plt.show()
