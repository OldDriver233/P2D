import numpy as np
import matplotlib.pyplot as plt

y = np.loadtxt('build/out.csv', delimiter=',')

plt.plot(y[:, 0], y[:, 1])
plt.show()