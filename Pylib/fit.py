################################
#			       #
#  Library for Nice Fits       #
#			       #
################################

import numpy as np
from scipy.optimize import curve_fit

#def chi2 functions
def chi2(y_measure,y_predict,errors):
    return np.sum( (y_measure - y_predict)**2 / errors**2 )

def chi2reduced(y_measure, y_predict, errors, number_of_parameters):
    return chi2(y_measure, y_predict, errors)/(y_measure.size - number_of_parameters)

#fit procedure 
def FitData(func, x, y, yerr, *param):
    popt, pcov = curve_fit(func, x, y, p0 = param, sigma = yerr, absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    chi2r = chi2reduced(y, func(x, *popt), yerr, len(param))
    
    #print the results
    for i in range(len(popt)):
        print(str(i+1) + str("_param"), " = ", popt[i], "+-", perr[i])
    
    print("chi2reduced = ", chi2r, "    ddof = ", len(y)-len(popt))
    print("----------------------------------------------------")
    
    return popt, perr, chi2r

#print functions
def PrintParams(params, perr, chi2red, num_data, savefig=False, namefile="fit_output"):
    if savefig:
        for index in range(len(params)):
            print(params[index], perr[index], end = " ", file=open(namefile, "a"))
        print(chi2red, num_data-len(params), file=open(namefile, "a"))
    else:
        for index in range(len(params)):
            print(params[index], perr[index], end = " ")
        print(chi2red, num_data-len(params))

