import matplotlib
matplotlib.use("Qt4Agg")

from plotqd import reader

import sys
from scipy.linalg import lstsq 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# read in data
data = reader(sys.argv[1])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

points = ax.scatter(xs=data[:,4], ys=data[:,5], zs=data[:,1], s=data[:,1], cmap=plt.cm.viridis)
    
deg = int(sys.argv[2]) + 1

def f(xi, i):
    """ function form of each LS term """
    return xi**i
# end function fcos

def xFun(x, *n):
    """ general fit function """
    X = np.zeros(np.shape(x))
    for i in xrange(len(X)):
        for j in xrange(len(n)):
            X[i] += n[j] * f(x[i], j)
        # end forj
    # end fori
    return X
# end function xFun

# fit (potentially non-linear) polynomial
Xmat = np.c_[np.ones(np.shape(data[:,4:6])[0]), data[:,4:6],
        np.prod(data[:,4:6], axis=1), data[:,4:6]**2]
C,_,_,_ = lstsq(Xmat, data[:,1]+data[:,3])


X, Y = np.meshgrid(data[:,4], data[:,5])
XX = X.flatten()
YY = Y.flatten()
Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2],
        C).reshape(X.shape)

ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.5)

plt.show()
