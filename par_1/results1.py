import numpy as np
import matplotlib.pyplot as plt

plt.plot([1, 2, 4, 5, 10, 20, 30, 40], 
         [14.5236, 23.0235, 50.1851, 56.7034, 107.885, 169.457, 217.823, 246.483])

plt.plot([1, 2, 4, 5, 10, 20, 30, 40], 
         14.5236 * np.array([1, 2, 4, 5, 10, 20, 30, 40]))

plt.show()