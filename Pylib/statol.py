################################
#			       #
#  Library for Data Analysis   #
#			       #
################################

import numpy as np



# takes in input np.array(data) annd a size N. It returns the blocked array with block length N
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


def block_sum(data, N):
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
	
    return np.sum(shaped_data, axis = 1)


def jack_for_single(func, data, N):
    if N <= 0:
        print("ERROR: the size of the block must be >0!")
        exit(1)
        
    num_blocks = int(len(data)/N)

    if num_blocks <=0:
        print("ERROR: the blocksize is too big!")
        exit(1)
        
    to_be_cut = len(data)%N
    cutted_data = data[to_be_cut:]
    
    eval_data = func(cutted_data)
 
    mean_value = np.mean(eval_data)
    
    # jacknife error procedure
    sums = block_sum(eval_data, N)
    total_sum = np.sum(sums)
    repeted_total_sum = np.tile(total_sum, num_blocks)
    jack_samples = (repeted_total_sum - sums)/((num_blocks - 1)*N)
    
    return mean_value, np.std(jack_samples)*np.sqrt(num_blocks - 1)
    

def jack_for_composite(comp_func, N, *args):
    """ Compute the jacknife analysis on comp_func(<func(data)>) 
        args is a list composed by list of the form [func, data].
        In the analysis N is used as the block size"""
    
    if N <= 0:
        print("ERROR: the size of the block must be >0!")
        exit(1)
        
    # jakcnife error list
    jack_list = []
    
    # mean values list
    mean_values_list = []
    
    for el in args:
        func, data = el
        
        num_blocks = int(len(data)/N)

        if num_blocks <=0:
            print("ERROR: the blocksize is too big!")
            exit(1)
        
        to_be_cut = len(data)%N
        cutted_data = data[to_be_cut:]
    
        eval_data = func(cutted_data)
 
        mean = np.mean(eval_data)
    
        # jacknife error procedure
        sums = block_sum(eval_data, N)
        total_sum = np.sum(sums)
        repeted_total_sum = np.tile(total_sum, num_blocks)
        jack_samples = (repeted_total_sum - sums)/((num_blocks - 1)*N)
        
        jack_list.append(jack_samples)
        mean_values_list.append(mean)        

    composite_jack = comp_func(jack_list)
    composite_jack = np.array(composite_jack)
    
    mean_value = num_blocks*comp_func(mean_values_list) - (num_blocks - 1)*np.mean(composite_jack)   
    err = np.std(composite_jack, ddof=1)*np.sqrt(num_blocks - 1)
    
    return mean_value, err
    


    
    
    
        
    
