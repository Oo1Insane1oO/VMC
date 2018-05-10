import numpy as np
import matplotlib.pyplot as plt
import sys

data = np.loadtxt(sys.argv[1])

r = np.linspace(0,3,500)

plt.plot(r, data, 'o', markersize=0.5)
plt.show()
