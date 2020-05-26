###############################
#							  #
#  Library for Data Analysis  #
#							  #
###############################
import numpy as np



# takes in input np.array(data) e returns the blocked array with block length N
def block_data(data, N):
	if N <= 0:
		print("ERROR: the size of the block must be > 0!")
		exit(1)
    
	blocked_lenght = int(len(data)/N)
    
	if blocked_lenght <= 0:
		print("ERROR: the  blocksize is too big!")
		exit(1)

	to_be_cut = len(data)%N
	cutted_data = data[to_be_cut:]
	shaped_data = cutted_data.reshape(blocked_lenght, N)
	
	return np.mean(shaped_data, axis = 1)
