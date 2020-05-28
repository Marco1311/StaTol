import numpy as np
import statol as stat
from random import * 


therm = 10000
delta = 0.1

#### Gaussian calculator ####
mu = 5
sigma = 1

def gauss(x):
	return np.exp(-(x-mu)**2/(2*(sigma**2)))

#### Metropolis ####
sample=[]
sample.append(0)
for i in range(10000000):
	if sample[i]-delta >0:
		candidate = sample[i] - delta + random()*2*delta
	else:
		candidate = sample[i] + random()*delta
	P_old = gauss(sample[i])	
	P_new = gauss(candidate)

	coin = random()
	if P_new > P_old:
		sample.append(candidate)
	else:
		if  coin < P_new/P_old:
			sample.append(candidate)
		else:
			sample.append(sample[i])
#	if i > therm:
#		print(sample[i+1])  


sample = np.array(sample)

for i in range(15):
	aux=stat.block_data(sample, 2**(i))
	print(np.mean(aux), np.std(aux, ddof=1)/np.sqrt(len(aux)), len(aux))
