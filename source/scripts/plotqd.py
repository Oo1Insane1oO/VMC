import matplotlib
matplotlib.use("Qt4Agg")

import numpy as np
from scipy.optimize import curve_fit 
import matplotlib.pyplot as plt
import os
import sys
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std

np.set_printoptions(linewidth=1000)

def arrayAllEqual(array):
    """ returne Ture if all values in array are equal to eachother """
    t = False
    for i in xrange(1,len(array)):
        if array[i] != array[i-1]:
            t = True
            break
        # end if
    # end fori
    return t
# end function arrayAllEqual

def reader(directory):
    data = []
    for root, dirs, files in os.walk(directory):
        """ grab root directory, directory and files """
        for f in files:
            """ loop over files  """
            if f[-4:] == ".txt":
                with open(os.path.join(root, f), 'r') as openfile:
                    """ open file for reading """
                    tmp = [] # temporary list for values
                    for line in openfile:
                        """ loop over lines in file and split """
                        words = line.split()
                        if len(words) != 0:
                            """ only append value from non-empty lines """
                            tmp.append(float(words[1]))
                        # end if
                    # end for
                    data.append(tmp) # append filled tmp into data list
                # end with open
            # end if
        # end for
    # and for

    # convert to numpy array and return sorted in increasing order with respect
    # to parameters
    data = np.array(data);
    sortIndices = []
    for i in xrange(4, np.shape(data)[1]):
        """ loop over parameters and find which columns need to be sorted """
        if arrayAllEqual(data[:,i]):
            """ only take columns were values are not equal to eachother """
            sortIndices.append(i)
        # end if
    # end fori

    return data[np.lexsort([data[:,i] for i in sortIndices])]
# end function

def svdsolve(A, rhs):
    """ solve linear system Ax = rhs for x with svd decomposition """

    # make sure matrix and rhs is aligned
    assert np.shape(A)[0] == np.shape(rhs)[0], ("A and rhs not aligned,"
            "shape(A): %s, shape(rhs): %s" % (str(np.shape(A)),
                str(np.shape(rhs))))

    # solve using svd decomposition and take care of extensions in case matrix
    # in non-square
    u, s, v = np.linalg.svd(A, full_matrices=(True if
        (np.shape(A)[0]==np.shape(A)[1]) else False))
    return np.dot(v, np.divide(np.dot(u.transpose(), rhs), s))
# end function svdsolve

def setLSmatrix(x, n, f):
    """ set matrix in LS problem """
    n += 1
    A = np.zeros((len(x), n))
    for i in xrange(len(x)):
        for j in xrange(n):
            A[i,j] = f(x[i], j)
        # end forj
    # end fori
    return A
# end function setLSmatrix

if __name__ == "__main__":
    """ standard 1-parameter run """

    # read in data and assign to array
    data = reader(sys.argv[1])

    # plot with errorbar
    plt.errorbar(data[:,4], data[:,1], yerr=data[:,3], fmt='o', markersize=1.5,
            ecolor='blue', capthick=1, label="Data-points", alpha=0.9)

    # fill area between errorbars
    plt.fill_between(data[:,4], data[:,1]-data[:,3], data[:,1]+data[:,3],
            alpha=0.2, facecolor='blue')

    # set title and labels
    plt.title("$\omega = %s$" % str(data[0,0]))
    plt.xlabel("$\\alpha$")
    plt.ylabel("$E$")

    ### perform regression analysis ###

#     # create matrix (predictor stacked with response)
#     X = sm.add_constant(np.column_stack((data[:,4], data[:,1])), prepend=False)
# 
#     # calculate fit with ordinary least squares
#     estimated = sm.OLS(data[:,1], X).fit()
# 
#     # find score and error in least square
#     prstd, errl, erru = wls_prediction_std(estimated)
# 
#     # plot prediction
#     plt.errorbar(data[:,4], estimated.predict(X), yerr=prstd, fmt='-',
#             ecolor='red', color='red', label="OLS fit", alpha=0.5)
#     plt.legend(loc='best')

    deg = int(sys.argv[2]) + 1

    def f(xi, i):
        """ function form of each LS term """
        return xi**i
    # end function fcos

    def xFun(x, *n):
        """ general fit function """
        X = np.zeros(np.shape(x))
        for i in xrange(x.size):
            for j in xrange(len(n)):
                X[i] += n[j] * f(x[i], j)
            # end forj
        # end fori
        return X
    # end function xFun

    # fit (potentially non-linear) polynomial
    P, C = curve_fit(xFun, data[:,4], data[:,1] + data[:,3], p0=tuple([0 for i
        in range(deg)])) 

    # make plot
    plt.errorbar(data[:,4], xFun(data[:,4], *P), yerr=0,
            capthick=1, markersize=1.5, ecolor='red', color='red', fmt='-',
            label=("OLS fit, deg=%i" % deg), alpha=0.5)
    plt.legend(loc='best')

    plt.show()
# end if main
