import numpy as np
import statol as stat
from random import random 


therm = 10000
delta = 0.1

#### Gaussian calculator ####
mu = 0
sigma = 1

def gauss(x):
	return np.exp(-(x-mu)**2/(2*(sigma**2)))

#### Metropolis ####
sample=[]
sample.append(0)
for i in range(10000000):
    candidate = sample[i] - delta + random()*2*delta
    
    P_old = gauss(sample[i])	
    P_new = gauss(candidate)

    coin = random()
    if P_new > P_old:
        sample.append(candidate)
    else:
       if coin < P_new/P_old:
           sample.append(candidate)
       else:
           sample.append(sample[i])
	#f i > therm:
	#	print(sample[i+1], file=open('MC_story.dat', "a"))  


sample = np.array(sample)
sample = sample[10*therm:]

def x4(data):
    return data*data*data*data

def x2(data):
    return data*data

def binder(args):
    return args[0]/(3*(args[1]**2))

for i in range(20):
    mean, err = stat.jack_for_composite(binder, 2**i, [x4, sample], [x2, sample])
    print(mean, err)
    #aux=stat.block_data(sample, 2**(i))
	#print(np.mean(aux), np.std(aux, ddof=1)/np.sqrt(len(aux)), len(aux))