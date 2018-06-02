import numpy as np
import matplotlib.pyplot as plt
import sys

fname = sys.argv[1]

data = np.loadtxt(fname)

rmaxIdx = fname.index("x")
rmax = float(fname[rmaxIdx+1:rmaxIdx+4])

r = np.linspace(0,rmax,600)

plt.plot(r, data, 'o', markersize=0.5)
plt.show()
